
MODULE SEMILAG

USE NUMERICS 
USE UTILS 


IMPLICIT NONE

REAL(RP),ALLOCATABLE, DIMENSION(:,:,:) :: U_SL
REAL(RP),ALLOCATABLE, DIMENSION(:,:) :: E_SL,ENERGY_SL


CONTAINS

SUBROUTINE SL_FIELD(FUN,A) 
  
  REAL(RP),DIMENSION(:,:),INTENT(IN) :: FUN
  REAL(RP),DIMENSION(:),INTENT(OUT) , ALLOCATABLE :: A
  
  INTEGER ::K,I
  REAL(RP)::S

  ALLOCATE(A(NX))
  
  A(1)=ZERO
    DO I=1,NX-1

      A(I+1)=A(I)+(UN-QUAD(DV,FUN(I,:)))*DX
            
    END DO
   
   S=SUM(A)/NX
     
     DO I=1,NX-1
     A(I)=A(I)-S
        
     END DO
     
     A(NX)=A(1)
 RETURN
 
END SUBROUTINE  

SUBROUTINE ENERGY(N,TIM)

INTEGER, INTENT(IN) ::N
REAL(RP) , INTENT(IN) ::TIM
REAL(RP), DIMENSION(NX) ::TMPBIS
INTEGER ::I
REAL(RP) ::S


DO I=1,NX
    TMPBIS(I)=E_SL(N,I)**2
END DO

ENERGY_SL(N,1)=TIM
ENERGY_SL(N,2)=SQRT(QUAD(DX,TMPBIS)*0.5_RP)
RETURN
END SUBROUTINE ENERGY









SUBROUTINE SL_SOLVER


INTEGER :: I,J,P,K,M,MM
REAL(RP):: ALPHA,CHRONO

REAL(rp),ALLOCATABLE, DIMENSION(:) :: TMP
REAL(rp),ALLOCATABLE, DIMENSION(:,:) :: Q,U,U_STAR

ALLOCATE(Q(NX,NV))
ALLOCATE(U(NX,NV))
ALLOCATE(U_STAR(NX,NV))
ALLOCATE(E_SL(NT,NX))
ALLOCATE(U_SL(NT,NX,NV))
ALLOCATE(ENERGY_SL(NT,2))


DT=TMAX/NT

U_SL(1,:,:)=FINIT
CALL SL_FIELD(FINIT,TMP)
E_SL(1,:)=TMP

CHRONO=ZERO



WRITE(*,*)"CFL_X",ABS(DT*VMAX/DX),"CFL_V",ABS(DT*MAXVAL(TMP)/DV)
WRITE(*,*)'LE NOMBRE D ITERATIONS EST ',NT


DO K=1,NT-1 ! BOUCLE EN TEMPS 

!------------------------TRANSPORT EN VITESSE---------------------------              
      
      DO I=1,NX-1
              
              P=FLOOR(-E_SL(K,I)*DT*0.5_RP/DV ) 
              ALPHA=ABS(-E_SL(K,I)*DT*0.5_RP/DV- P)
             
              DO J=2,NV-1 
               
                IF(J-P <= 1 .OR. J-P > NV )THEN
                   U(I,J)=ZERO
                    
                ELSE   
                   U(I,J)=(UN-ALPHA)*U_SL(K,I,J-P)+ALPHA*U_SL(K,I,J-P-1)   
                END IF
                   
              END DO
               
        END DO
      
        U(:,1)=ZERO
        U(:,NV)=ZERO
        U(NX,:)=U(1,:)
       
        

!---------------------TRANSPORT EN ESPACE-------------------------------      
      
       DO J=2,NV-1
       
              P=FLOOR(V(J)*DT/DX)
              ALPHA=ABS(V(J)*DT/DX - P  )

              DO I=1,NX-1
              
                
                IF(I-P==0 .OR. I-P==1 )THEN
               
                Q(I,J)=U(I,J)
                ELSE
               
                Q(I,J)=(UN-ALPHA)*U(I-P,J)+ ALPHA*U(I-P-1,J)
                END IF                                
              END DO
              
        END DO
        
        Q(:,1)=ZERO
        Q(:,NV)=ZERO
        Q(NX,:)=Q(1,:)


!-----------------------CALCUL DU CHAMP E-------------------------------        
        
        !CALL SL_FIELD(Q,TMP)
        
        
        
!------------------------TRANSPORT EN VITESSE---------------------------              
      
      DO I=1,NX-1
        
              P=FLOOR(-TMP(I)*DT*0.5_RP/DV)
              ALPHA=ABS(-TMP(I)*DT*0.5_RP/DV- P)
                        
              DO J=2,NV-1 
                  
					IF(J-P <= 1 .OR. J-P > NV )THEN
					
	                    U_SL(K,I,J)=ZERO
	                  
	                   
	                ELSE
	                
						U_SL(K+1,I,J)=(UN-ALPHA)*Q(I,J-P)+ ALPHA*Q(I,J-P-1)
	               
                    END IF  
 
                 
              END DO
              
        END DO

        U_SL(K+1,:,1)=ZERO
        U_SL(K+1,:,NV)=ZERO
        U_SL(K+1,NX,:)=U_SL(K+1,1,:)


!----------------------ACTUALISATION DE E-------------------------------
      
        CALL SL_FIELD(U_SL(K+1,:,:),TMP)
        
        E_SL(K+1,:)=TMP
        CHRONO=CHRONO+DT
        CALL ENERGY(K+1,CHRONO)

END DO
 
CALL SAVE_RESULT(NT)

WRITE(*,*) "LE TEMPS FINAL EST ", CHRONO

CALL SAVE_SL_ENERGY
CALL SPACE_ERR(U_SL(2,:,:),DT)

 
DEALLOCATE(Q,U,U_STAR,E_SL,U_SL,ENERGY_SL)
        

END SUBROUTINE SL_SOLVER

SUBROUTINE SAVE_RESULT(K)

 
   
 INTEGER , INTENT(IN) ::K
 integer :: I,J
 !stockage dans un fichier
    
      open(unit=14, file='SL_LIN.txt', status='REPLACE')
        
       DO I=1,SIZE(FINIT,1)
            DO J=1,SIZE(FINIT,2)
    
                    WRITE(14,*)  X(I),' ',V(J),' ',U_SL(K,I,J)

    
            END DO     
       END DO
    
       close(14)
   
END SUBROUTINE save_result



SUBROUTINE SAVE_SL_ENERGY

 
   
INTEGER :: I,K
REAL(RP), DIMENSION(NX) ::TMP
REAL(RP) ::S,TIM

TIM=ZERO
DO K=1,NT

    DO I=1,NX
        TMP(I)=E_SL(K,I)**2
    END DO

    ENERGY_SL(K,1)=TIM
    ENERGY_SL(K,2)=SQRT(QUAD(DX,TMP)/2._RP)
    
    TIM=TIM+DT
END DO 
    
  open(unit=17, file='SL_ENERGY.txt', status='REPLACE')
        
       DO I=1,NT

                    WRITE(17,*)  ENERGY_SL(I,1),'  ', LOG(ENERGY_SL(I,2))
 
       END DO
    
       close(17)
   
  END SUBROUTINE SAVE_SL_ENERGY





END MODULE SEMILAG
