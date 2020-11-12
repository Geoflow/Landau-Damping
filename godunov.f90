MODULE GODUNOV
USE NUMERICS 
USE UTILS 
IMPLICIT NONE 

REAL(RP),ALLOCATABLE, DIMENSION(:,:,:) :: U_GS
REAL(RP),ALLOCATABLE, DIMENSION(:,:) :: E_GS,ENERGY


CONTAINS

SUBROUTINE ENERGY_FD(N,TIM)

INTEGER, INTENT(IN) ::N
REAL(RP) , INTENT(IN) ::TIM
REAL(RP), DIMENSION(NX) ::TMPBIS
INTEGER ::I
REAL(RP) ::S


DO I=1,NX
    TMPBIS(I)=E_GS(N,I)**2
END DO

ENERGY(N,1)=TIM
ENERGY(N,2)=SQRT(QUAD(DX,TMPBIS)*0.5_RP)
RETURN
END SUBROUTINE ENERGY_FD

SUBROUTINE E_FIELD(N,A) 
  
  INTEGER, INTENT(IN) ::N
  REAL(RP),DIMENSION(:),INTENT(OUT) , ALLOCATABLE :: A
  
  INTEGER ::I
  REAL(RP)::S

  
  ALLOCATE(A(NX))
  
  A(:)=zero
    DO I=2,NX-1

      A(I+1)=A(I)+(UN-QUAD(DV,U_GS(N,I,:)))*DX
            
    END DO
     
     S=SUM(A)/NX
     
     DO I=1,NX
     A(I)=A(I)-S
        
     END DO
     
     A(NX)=A(1)
 RETURN
END SUBROUTINE  


SUBROUTINE GS    ! GODUNOV SPLITTING SCHEMA UPSTREAM

REAL(rp),ALLOCATABLE, DIMENSION(:,:) :: Q,U,U_STAR


REAL(rp),ALLOCATABLE, DIMENSION(:) :: TMP
INTEGER ::I,J,K
REAL(RP)::CHRONO

ALLOCATE(Q(NX,NV))
ALLOCATE(U(NX,NV))
ALLOCATE(TMP(NX) )
ALLOCATE(ENERGY(NT,2))
ALLOCATE(E_GS(NT,NX))
ALLOCATE(U_GS(NT,NX,NV))

DT=TMAX/NT

U_GS(1,:,:)=FINIT
Q=FINIT
CALL E_FIELD(1,TMP)

E_GS(1,:)=TMP
     
WRITE(*,*)"CFL_X",ABS(DT*10._rp/DX),"CFL_V",ABS(DT*maxval(tmp)/DV)
WRITE(*,*)'LE NOMBRE D ITERATIONS EST ',NT


CHRONO=ZERO

!---------------------------BOUCLE EN TEMPS-------------------------

DO K=1,NT-1


!--------------------------TRANSPORT EN VITESSE-------------------------      
        
        DO J=2,NV-1
				DO I=1,NX-1
                 
                  U(I,J)=Q(I,J)-(MAX(ZERO,-TMP(I))*( Q(I,J)-Q(I,J-1))+MIN(ZERO,-TMP(I))*(Q(I,J+1)-Q(I,J)))*DT/DV
                   
                END DO
        END DO
        
        U(NX,:)=U(1,:)
        
   
!--------------------------TRANSPORT EN ESPACE ------------------------       
        U(:,1)=ZERO
        U(:,NV)=ZERO
        
        DO J=2,NV-1
        
            Q(1,J)=U(1,J)-V(J)*DT*(U(2,J)-U(1,J))/DX
            
            DO I=2,NX-1
            
                Q(I,J)=U(I,J)-(MAX(ZERO,V(J))*(U(I,J)-U(I-1,J))/DX+MIN(ZERO,V(J))*(U(I+1,J)-U(I,J))/DX)*DT
                                                
            END DO
        END DO
      
        Q(NX,:)=Q(1,:)
        Q(:,1)=ZERO
        Q(:,NV)=ZERO
        
!---------------------------CALCUL DU CHAMP E---------------------------
       
        U_GS(K+1,:,:)=Q     
        CALL E_FIELD(K+1,TMP)
        
        E_GS(K+1,:)=TMP 
        CALL ENERGY_FD(K+1,CHRONO)
        CHRONO=CHRONO+DT
      
END DO

!-----------------------------SAUVEGARDE---------------------------------

        CALL SAVE_RESULT(NT)
        CALL SAVE_ENERGY(NT)
        
        WRITE(*,*) "LE TEMPS FINAL EST ", CHRONO
        CALL SPACE_ERR(U_GS(2,:,:),DT)

        DEALLOCATE(U,Q,ENERGY,U_GS,TMP)
        
END SUBROUTINE GS

SUBROUTINE SAVE_RESULT(K)

 
   
 INTEGER , INTENT(IN) ::K
 integer :: I,J
 !stockage dans un fichier
    
      open(unit=13, file='godunov_split.txt', status='REPLACE')
        
       DO I=1,SIZE(FINIT,1)
            DO J=1,SIZE(FINIT,2)
    
                    WRITE(13,*)  X(I),' ',V(J),' ',U_GS(K,I,J)

    
            END DO     
       END DO
    
       close(13)
   
END SUBROUTINE save_result

SUBROUTINE SAVE_ENERGY(N)

 
   
 integer , INTENT(IN):: N
 integer :: I
 !stockage dans un fichier
    
      open(unit=17, file='GD_ENERGY.txt', status='REPLACE')
        
       DO I=2,NT

                    WRITE(17,*)  ENERGY(I,1),' ', LOG(ENERGY(I,2))
 
       END DO
    
       close(17)
   
  END SUBROUTINE SAVE_ENERGY


END MODULE GODUNOV
