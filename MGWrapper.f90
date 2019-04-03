! compile with "ifort -c MGWrapper.f90 "
! create a library with "ar rcs MGWrapper.a MGWrapper.o /users/pep/mschulze/lib/MadGraph-2.49/switchmom.o /users/pep/mschulze/lib/HELAS-3.0/*.o"

! the Madgraph file has to be compiled with "ifort -c -I/mnt/pep/mschulze/lib/HELAS-3.0 emep_taptam.f" 
! 
! the HELAS library has to be compiled with "ifort -lifcore -limf ....."

function MGWrapper(N,Mom,Res) BIND(C)
use, intrinsic :: ISO_C_BINDING 
implicit none
integer(C_INT) :: N
real(C_DOUBLE) :: Mom(1:N,1:4),Res,MGWrapper


      call coupsm(0)
      call SEMEP_TAPTAM(Mom,Res)
      

MGWrapper=1.0
return 
end function 

! 
! call this function from C with 
!    int NPart=4;
!    double MGRes;
!    double MomExtMG[4][4];    \\ first: particle index; second: Lorentz index
!    mgwrapper(&NPart,MomExtMG,&MGRes);
!    printf("MG Result %f \n",MGRes); 
