program TUCKEL

!_______________________________________________________________________________
!
!                       Welcome to TUCKEL
!                  Another Hückel implementation
!                               by
!                   Bustamante Carlos Mauricio
!                     carlosmbqca@gmail.com
!
! TUCKEL is a new implementation of Extended Hückel Molecular Orbital theory
! (EHMO for short), which use is destined to university teaching.
!
! Copyright (C) 2015  Bustamante Carlos Mauricio
!
! This file is part of TUCKEL.
!
! TUCKEL is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! TUCKEL is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! A copy of the licence can be found in the root directory of
! TUCKEL.  If not, see <http://www.gnu.org/licenses/>.
!
!_______________________________________________________________________________

   implicit none

!****Vatiables:

   integer:: i, j, k, l                                               !counters
   integer:: e, nb, t                                                 !flags
   integer:: nro                                                      !number of atoms
   integer:: nf                                                       !number of basis
   integer:: ele                                                      !flag for electrons
   integer:: electrons                                                !number of electrons
   integer:: q                                                        !charge
   integer:: option1, option2                                         !options of menu
   integer:: nro_shells                                               !number of shells
   integer:: orbitals                                                 !number of orbitals per atom
   integer:: LWORK, INFO                                              !lapack subroutine arguments
   integer:: step_min, step_max, steps, atom_ref                      !Arguments needed for scan process
   integer:: ios                                                      !Value that indicate the state of open files
   character (len=1) :: atom                                          !atom's name flag
   character (len=1) :: C="C", H="H", B="B", O="O", N="N", F="F"      !arom's name
   character (len=20) :: indoc                                        !name of the input
   character (len=30) :: doc                                          !name + extension
   character (len=4) :: extin= ".dat", extout=".out"                  !extension
   character (len=100) :: filename                                    !filename
   character (len=100) :: coordenada                                  !kind of coordinate of scan
   character (len=20) :: str                                          !name of function
   real*8:: x1, y1, z1, x2, y2, z2, a1, a2, d1, d2                    !flags for coordinates, alpha and d
   real*8:: s2, sp, p2, pp                                            !names of functions
   real*8:: total                                                     !flag for overlap elements
   real*8:: normv                                                     !normalization value for eigenvectors
   real*8:: band                                                      !flag for elemnts of P matrix
   real*8:: netq                                                      !flag for P diagonal elements
   real*8:: totalenergy                                               !total energy
   real*8:: delta                                                     !scan parameter
   real*8:: inicio                                                    !initial value of scan
   real*8, dimension(3):: ri                                          !read array of coordinates
   real*8:: a=0.5291772089                                            !Bohr radius
   real*8:: hartree=0.036749309                                       !convertion value eV -> Hartree


!*************************************Interfaces to indicate the dynamics arrays************************************

   interface
      subroutine element (atom, gi, ai, di, orbitals)
         character (len=1), intent(in) :: atom
         character (len=4) :: C="C", H="H", B="B", O="O", N="N", F="F"
         character(len=4), allocatable, intent(out):: gi(:, :)
         real*8, allocatable, intent(out):: ai(:, :), di(:, :)
         integer, intent(out):: orbitals
      end subroutine element

      subroutine z_elements(distance, simb, ang, die, union, nro)
         real*8, allocatable, intent(inout):: distance(:), ang(:), die(:)
         integer, allocatable, intent(inout):: union(:, :)
         character (len=4), allocatable, intent(inout):: simb(:)
         integer, intent(in):: nro
         integer:: i
      end subroutine z_elements

      subroutine z_translate(distance, ang, die, union, r, nro)
         real*8, allocatable, intent(in):: distance(:), ang(:), die(:)
         real*8, allocatable, intent(inout):: r(:, :)
         real*8:: v1, v2, v3, r1(3), r2(3), r3(3), x, y, z, nule=0
         integer, allocatable, intent(in):: union(:, :)
         integer, intent(in):: nro
         integer:: i
      end subroutine z_translate

     subroutine S_H_matrix (S, ham, simb, g, alfa, d, r, iden, norm, nf)
        real*8, allocatable, intent(inout):: S(:, :), norm(:), ham(:, :)
        real*8, allocatable, intent(in):: d(:, :), alfa(:, :), r(:, :)
        character (len=4), allocatable, intent(in):: g(:, :), simb(:)
        integer, allocatable, intent(in):: iden(:)
        integer, intent(in):: nf
        integer:: i, j, k, l
        real*8:: s2, sp, p2, pp
        real*8:: x1, y1, z1, a1, a2, d1, d2, total
     end subroutine S_H_matrix

      
   end interface
      
!********************************************************************End of interface

!***Allocatable variables

   character(len=4), allocatable:: gi(:, :)                           !flag of g array
   character(len=4), allocatable:: g(:, :)                            !contains the identity of each kind gaussian primitive
   character(len=4), allocatable:: simb(:)                            !array of atoms simbols
   real*8, allocatable:: ai(:, :), di(:, :)                           !flags of alfa and d array
   real*8, allocatable:: alfa(:, :)                                   !array that contains the primitive coefficients
   real*8, allocatable:: d(:, :)                                      !array that contains the contraction coefficients
   real*8, allocatable:: r(:, :)                                      !array that contains the cartesian coordinates
   real*8, allocatable:: S(:, :)                                      !Overlap matrix
   real*8, allocatable:: ham(:, :)                                    !Hamiltonian matrix
   real*8, allocatable:: m1(:, :)                                     !array that will contain the eigenvectors after diagonalisation
   real*8, allocatable:: m2(:, :)                                     !array that will contain S matrix previous diagonalisation
   real*8, allocatable:: energy(:)                                    !array that contain the energy of each state
   real*8, allocatable:: norm(:)                                      !narmalization array for S matrix
   real*8, allocatable:: P(:, :)                                      !density array
   real*8, allocatable:: orden(:, :)                                  !array that contains bonds orders and electron densities
   real*8, allocatable:: nval(:)                                      !array that contains valnces of each atom
   real*8, allocatable:: WORK(:)                                      !lapack subroutine argument
   real*8, allocatable:: energy_scan(:)                               !array that contain the energy of each scan
   real*8, allocatable:: distance(:), ang(:), die(:)                  !parameters of Z-matrix
   real*8, allocatable:: r_cc(:), prim_exp(:), contr_coeff(:), mo_coeff(:) !store arrays for write values in "fchk" file
   integer, allocatable:: iden(:)                                     !array that contain the atom number where belong a base
   integer, allocatable:: z_atom(:)                                   !array that contain the atomic number of each atom
   integer, allocatable:: shell_type(:), shell_iden(:)                !store arrays for write values in "fchk" file
   integer, allocatable:: union(:, :)                                 !array that contains the atom numbers assoiciated to a Z-coordenate

!***********************Menu********************************************************

print*, "***Bienvenido a TUCKEL***"
print*, "------------------------------------------------------------------"
print*, "TUCKEL  Copyright (C) 2015  Bustamante Carlos Mauricio"
print*, "------------------------------------------------------------------"
print*, "Por favor elija una de las opciones:"
print*, "1-Calculo de Energías, densidades electronicas y ordenes de enlace"
print*, "2-Scan de coordenadas internas"
print*, "3-Salir del programa"
read*, option1

if (option1==1) then
   print*, "Por favor elija el tipo de coordenadas del input:"
   print*, "1-Cartesianas"
   print*, "2-Matriz-Z"
   read*, option2
else if (option1==3) then
   stop
end if

!*********************Selection and opening of the input***************************
   print*, "Inserte el nombre del archivo sin la extensión .dat :"
   read*, indoc
   indoc=adjustr(indoc)
   doc= indoc // extin
   doc=adjustl(doc)
   open (unit=13, file=doc, status="old", action="read", iostat=ios)
   if (ios/=0) then
      print*, "Error al abrir el archivo:"
      print*, doc
      stop
   end if

!*********************Counting the number of atoms*********************************
nro=0
do
   read(unit=13, fmt="(a)", iostat=ios) atom
   if (ios<0) then
      exit
   else if (ios>0) then
      print*, "Error al leer el input"
      stop
   else
      if (index(atom, C)>0 .or. index(atom, H)>0 .or. index(atom, N)>0 .or. index(atom, O)>0 .or. index(atom, B)>0 &
      .or. index(atom, F)>0) then
         nro=nro+1
      end if

   end if
end do
close (unit=13)

!*********************Counting the number of basis, shells and electrons******************************
   open (unit=13, file=doc, status="old", action="read", iostat=ios)
   if (ios/=0) then
      print*, "Error al abrir el archivo:"
      print*, doc
      stop
   end if
   
   read(13, *) filename
   read(13, *) q
   nf=0
   electrons=0
   ele=0
   nb=0
   nro_shells=0
   do i=1, nro
   read (13, *) atom
     
      if (atom == C) then
         nro_shells=nro_shells+2
         nb=4
         ele=4
      else if (atom == H) then
         nro_shells=nro_shells+1
         nb=1
         ele=1
      else if (atom == B) then
         nro_shells=nro_shells+2
         nb=4
         ele=3
      else if (atom == N) then
         nro_shells=nro_shells+2
         nb=4
         ele=5
      else if (atom == O) then
         nro_shells=nro_shells+2
         nb=4
         ele=6
      else if (atom == F) then
         nro_shells=nro_shells+2
         nb=4
         ele=7
      end if
      ele=ele+electrons
      electrons=ele
      nb=nb+nf
      nf=nb
   end do
   electrons=electrons-q

   if (option1==2) then                                  !If we want to make a scan at the end of the file we have to indicate:
      read(13, *) coordenada                             !The coordinate that we want scan: Enlace, Angulo o Diedro
         coordenada=trim(coordenada)
      read(13, *) atom_ref, delta, step_min, step_max    !Which atom, deltas of this coordinates, min step, and max step of deltas
         steps=step_max-step_min+1
         allocate(energy_scan(steps))
   end if

   close(unit=13)

!******************Allocating arrays**************************************

allocate(g(6, nf), alfa(6, nf), d(6, nf), iden(nf), r(3, nro), ham(nf, nf), norm(nf), simb(nro), nval(nro))
allocate(union(nro, 4), distance(nro), ang(nro), die(nro))
allocate(z_atom(nro), r_cc(3*nro), shell_type(nro_shells), prim_exp(6*nro_shells), contr_coeff(nro_shells*6), mo_coeff(nf**2))
allocate(shell_iden(nro_shells))


!******************Ordering all the data in arrays, and assigning the diagonal elements of Hamiltonian********
   
   open (unit=10, file=doc, status="old", action="read", iostat=ios)
   
   if (ios/=0) then
      print*, "Error to open file:"
      print*, doc
      stop
   end if

   read(10, *) filename

   read(10, *) q

   e=0

   if (option1==1 .and. option2==1) then                         !*** If the input have cartesian coordinates these will get in into r
      do i=1, nro            
         read (10, *) simb(i), ri(1), ri(2), ri(3)
         r(1:3, i)=ri(1:3)/a
      end do
   else if((option1==2).or.(option1==1 .and. option2==2)) then   !*** If the input have Z-matrix these coordinates will be ordered 
                                                                 !*** and then translated to cartesian coordinates
      union(1:nro, 1:3)=0
      ang(1:nro)=0
      die(1:nro)=0
      call z_elements(distance, simb, ang, die, union, nro)
      call z_translate(distance, ang, die, union, r, nro)

      r(1:3, 1:nro)=r(1:3, 1:nro)/a

   end if


   nb=0
   do i=1, nro                                                   !***Each atom is evaluated and is assigned to a atomic number, valence,
                                                                 !***diagonal element of H, and basis parameters
         call element (simb(i), gi, ai, di, orbitals)
         
      if (simb(i) == C) then
            nval(i)=4
            z_atom(i)=6
            shell_type(nb+1)=0
            shell_type(nb+2)=1
            shell_iden(nb+1:nb+2)=i
            nb=nb+2
            ham(e+1, e+1)=-21.4
            ham(e+2:e+4, e+2:e+4)=-11.4
         do j= e+1, e+4
            g(1:6, j)= gi(1:6,(j-e))
            alfa(1:6, j)= ai(1:6,(j-e))
            d(1:6, j)= di(1:6,(j-e)) 
            iden(j)=i
         end do
         e=e+4
      else if (simb(i) == B) then
            nval(i)=3
            z_atom(i)=5
            shell_type(nb+1)=0
            shell_type(nb+2)=1
            shell_iden(nb+1:nb+2)=i
            nb=nb+2
            ham(e+1, e+1)=-15.2
            ham(e+2:e+4, e+2:e+4)=-8.5
         do j= e+1, e+4
            g(1:6, j)= gi(1:6,(j-e))
            alfa(1:6, j)= ai(1:6,(j-e))
            d(1:6, j)= di(1:6,(j-e))
            iden(j)= i
         end do
         e=e+4
      else if (simb(i) == N) then
            nval(i)=5
            z_atom(i)=7
            shell_type(nb+1)=0
            shell_type(nb+2)=1
            shell_iden(nb+1:nb+2)=i
            nb=nb+2
            ham(e+1, e+1)=-25.58
            ham(e+2:e+4, e+2:e+4)=-13.9
         do j= e+1, e+4
            g(1:6, j)= gi(1:6,(j-e))
            alfa(1:6, j)= ai(1:6,(j-e))
            d(1:6, j)= di(1:6,(j-e))
            iden(j)= i
         end do
         e=e+4
      else if (simb(i) == O) then
            nval(i)=6
            z_atom(i)=8
            shell_type(nb+1)=0
            shell_type(nb+2)=1
            shell_iden(nb+1:nb+2)=i
            nb=nb+2
            ham(e+1, e+1)=-32.38
            ham(e+2:e+4, e+2:e+4)=-15.85
         do j= e+1, e+4
            g(1:6, j)= gi(1:6,(j-e))
            alfa(1:6, j)= ai(1:6,(j-e))
            d(1:6, j)= di(1:6,(j-e))
            iden(j)= i
         end do
         e=e+4
      else if (simb(i) == F) then 
            nval(i)=7
            z_atom(i)=9
            shell_type(nb+1)=0
            shell_type(nb+2)=1
            shell_iden(nb+1:nb+2)=i
            nb=nb+2
            ham(e+1, e+1)=-40.20
            ham(e+2:e+4, e+2:e+4)=-18.66
         do j= e+1, e+4
            g(1:6, j)= gi(1:6,(j-e))
            alfa(1:6, j)= ai(1:6,(j-e))
            d(1:6, j)= di(1:6,(j-e))
            iden(j)= i
         end do
         e=e+4
      else if (simb(i) == H) then
            nval(i)=1
            z_atom(i)=1
            shell_type(nb+1)=0
            shell_iden(nb+1)=i
            nb=nb+1
            ham(e+1, e+1)=-13.6
            j= e+1
            g(1:6, j)= gi(1:6, (j-e))
            alfa(1:6, j)= ai(1:6, (j-e))
            d(1:6, j)= di(1:6, (j-e))
            iden(j)= i
         e=e+1
      end if

   end do


!In this part the program make two different processes depending the first option choose


!***************************************Option 1: Energy, densities and bond order anlisis*****************************


allocate(S(nf, nf))

if (option1==1) then                                                          !*** If we select option 1

call S_H_matrix (S, ham, simb, g, alfa, d, r, iden, norm, nf)

LWORK=-1                                                                      !*** Using lapack libraries

allocate(m1(nf, nf), m2(nf, nf), energy(nf), WORK(LWORK))
   m1(1:nf, 1:nf)=ham(1:nf, 1:nf)
   m2(1:nf, 1:nf)=S(1:nf, 1:nf)

call dsygv (1, 'V', 'L', nf, m1, nf, m2, nf, energy, WORK, LWORK, INFO)
LWORK=WORK(1)
deallocate (WORK)
allocate (WORK(LWORK))
call dsygv (1, 'V', 'L', nf, m1, nf, m2, nf, energy, WORK, LWORK, INFO)

   if (INFO/=0) then
      print*, "Error during calculation of eigenvectors"
      stop
   end if

!*************************************Normalization of eigenvectors:
   do i=1, nf
      normv=0
      do j=1, nf
         normv=normv + m1(j, i)**2
      end do
      do j=1, nf
         m1(j, i)=m1(j, i)/(normv**(0.5))
      end do
   end do

!************************************Creating P array:

electrons=electrons/2

allocate (P(nf, nf))

P(1:nf, 1:nf)=0

   do i=1, nf
      do j=1, nf
         do k=1, electrons
            band=2*m1(i, k)*m1(j, k)
            P(i, j)=P(i, j) + band
         end do
      end do
   end do


!***********************************Creating orden array:

allocate (orden(nro, nro))

orden(1:nro, 1:nro)=0
      
   do i=1, nf
         band=P(i, i)
         orden(iden(i), iden(i))=orden(iden(i), iden(i))+band
   end do
band=0
   do i=1, nf
      do j= 1, nf
         if (iden(i)/=iden(j))then
            band=P(i, j)**2
            orden(iden(i), iden(j))=orden(iden(i), iden(j))+band
         end if
      end do
   end do


!************************************Creating an output file:

doc=indoc // extout
doc=adjustl(doc)

open(unit=11, file=doc, status="new", action="write", iostat=ios)
      if (ios/=0) then
      print*, "Error to create:"
      print*, doc
      stop
   end if

write(11, "(A)") "Gracias por usar TUCKEL"
write(11, "(A)") "--------------------------------------------------------"
write(11, "(A)") "TUCKEL  Copyright (C) 2015  Bustamante Carlos Mauricio"
write(11, "(A)") "--------------------------------------------------------"
write(11, "(A)")
write(11, "(A)") "**************************************************************"
write(11, "(A)") filename
write(11, "(A)") "**************************************************************"
write(11, "(A)") "Energia de los niveles ocupados"

totalenergy=0
do i=1, electrons
    write(11, "(A, F8.4, A)") "Estate"//trim(str(i))//":", energy(i), "eV"
    totalenergy=totalenergy + 2*energy(i)
end do
write(11, "(A, F15.4, A)"), "Energía total:", totalenergy, "eV"
write(11, "(A)") "******************************************************************"
write(11, "(A)") "Analisis de las densidades electronicas de Mulliken:"
   netq=0
   do i=1, nro
      write(11, "(A, F8.4)") adjustr(simb(i))//"("//trim(str(i))//"):", (nval(i)-orden(i, i))
      netq=netq+(nval(i)-orden(i, i))
   end do
write(11, "(A)") "******************************************************************" 
write(11, "(A, F8.4)") "Carga total:", netq
write(11, "(A)") "******************************************************************"
write(11, "(A)") "Ordenes de enlace de Wiberg:"
    do i=1, nro-1
      do j=i+1, nro
       write(11, "(A, F8.4)") adjustr(simb(i))//"("//trim(str(i))//")-"//trim(simb(j))//"("//trim(str(j))//"):", orden(i, j)
      end do
    end do

!********************************Creating a fchk file (for Avogadro):

doc=indoc // ".fchk"
doc=adjustl(doc)

open(unit=13, file=doc, status="new", action="write", iostat=ios)
      if (ios/=0) then
      print*, "Error to create:"
      print*, doc
      stop
   end if

write(13, "(A)") "Gracias por usar TUCKEL"
write(13, "(A)") "--------------------------------------------------------"
write(13, "(A)") "TUCKEL  Copyright (C) 2015  Bustamante Carlos Mauricio"
write(13, "(A)") "--------------------------------------------------------"
write(13, "(A)")
write(13, "(A)") filename
write(13, "(A)") "Number of atoms                            I                "//trim(str(nro))
write(13, "(A)") "Charge                                     I                "//trim(str(q))
write(13, "(A)") "Multiplicity                               I                1"
write(13, "(A)") "Number of electrons                        I               "//trim(str(electrons*2))
write(13, "(A)") "Number of alpha electrons                  I                "//trim(str(electrons))
write(13, "(A)") "Number of beta electrons                   I                "//trim(str(electrons))
write(13, "(A)") "Number of basis functions                  I               "//trim(str(nf))
write(13, "(A)") "Number of independent functions            I               "//trim(str(nf))
write(13, "(A)") "Number of point charges in /Mol/           I                0"
write(13, "(A)") "Number of translation vectors              I                0"
write(13, "(A)") "Atomic numbers                             I   N=           "//trim(str(nro))

   e=0
   nb=nro

   do while (nb/=0)
      if (nb>=5) then
         write(13, "(5I12)") (z_atom(i), i=e+1, e+5 )
         e=e+5
         nb=nb-5
      else
         write(13, "(5I12)") (z_atom(i), i=e+1, e+nb)
         exit
      end if
   end do

write(13, "(A)") "Current cartesian coordinates              R   N=           "//trim(str(nro*3))
   e=0
   do i=1, nro
      r_cc(e+1:e+3)=r(1:3, i)
      e=e+3
   end do

   e=0
   nb=nro*3
   do while (nb/=0)
      if (nb>=5) then
         write(13, "(5ES16.8E2)") (r_cc(i), i=e+1, e+5 )
         e=e+5
         nb=nb-5
      else
         write(13, "(5ES16.8E2)") (r_cc(i), i=e+1, e+nb)
         exit
      end if
   end do

write(13, "(A)") "Number of contracted shells                I                "//trim(str(nro_shells))
write(13, "(A)") "Number of primitive shells                 I               "//trim(str(nro_shells*6))
write(13, "(A)") "Shell types                                I   N=           "//trim(str(nro_shells))

   e=0
   nb=nro_shells
   do while (nb/=0)
      if (nb>=5) then
         write(13, "(5I12)") (shell_type(i), i=e+1, e+5 )
         e=e+5
         nb=nb-5
      else
         write(13, "(5I12)") (shell_type(i), i=e+1, e+nb)
         exit
      end if
   end do

write(13, "(A)") "Number of primitives per shell             I   N=           "//trim(str(nro_shells))
   e=0
   nb=nro_shells
   do while (nb/=0)
      if (nb>=5) then
         write(13, "(5I12)") (6, i=e+1, e+5 )
         e=e+5
         nb=nb-5
      else
         write(13, "(5I12)") (6, i=e+1, e+nb)
         exit
      end if
   end do

write(13, "(A)") "Shell to atom map                          I   N=           "//trim(str(nro_shells))

   e=0
   nb=nro_shells
   do while (nb/=0)
      if (nb>=5) then
         write(13, "(5I12)") (shell_iden(i), i=e+1, e+5 )
         e=e+5
         nb=nb-5
      else
         write(13, "(5I12)") (shell_iden(i), i=e+1, e+nb)
         exit
      end if
   end do

write(13, "(A)") "Primitive exponents                        R   N=          "//trim(str(nro_shells*6))

   e=0
   do i=1, nf
      if (g(1, i)=="ss" .or. g(1, i)=="px") then
         prim_exp(e+1:e+6)=alfa(1:6, i)
         e=e+6
      end if
   end do

   e=0
   nb=nro_shells*6
   do while (nb/=0)
      if (nb>=5) then
         write(13, "(5ES16.8E2)") (prim_exp(i), i=e+1, e+5 )
         e=e+5
         nb=nb-5
      else
         write(13, "(5ES16.8E2)") (prim_exp(i), i=e+1, e+nb)
         exit
      end if
   end do

write(13, "(A)") "Contraction coefficients                   R   N=          "//trim(str(nro_shells*6))

   e=0
   do i=1, nf
      if (g(1, i)=="ss".or. g(1, i)=="px") then
         contr_coeff(e+1:e+6)=d(1:6, i)
         e=e+6
      end if
   end do

   e=0
   nb=nro_shells*6
   do while (nb/=0)
      if (nb>=5) then
         write(13, "(5ES16.8E2)") (contr_coeff(i), i=e+1, e+5 )
         e=e+5
         nb=nb-5
      else
         write(13, "(5ES16.8E2)") (contr_coeff(i), i=e+1, e+nb)
         exit
      end if
   end do

write(13, "(A)") "Alpha Orbital Energies                     R   N=          "//trim(str(nf))

   e=0
   nb=nf
   do while (nb/=0)
      if (nb>=5) then
         write(13, "(5ES16.8E2)") (energy(i)*hartree, i=e+1, e+5 )
         e=e+5
         nb=nb-5
      else
         write(13, "(5ES16.8E2)") (energy(i)*hartree, i=e+1, e+nb)
         exit
      end if
   end do

write(13, "(A)") "Alpha MO coefficients                      R   N=          "//trim(str(nf**2))

   e=0
   do i=1, nf
      mo_coeff(e+1:e+nf)=m1(1:nf, i)
      e=e+nf
   end do

   e=0
   nb=nf**2
   do while (nb/=0)
      if (nb>=5) then
         write(13, "(5ES16.8E2)") (mo_coeff(i), i=e+1, e+5 )
         e=e+5
         nb=nb-5
      else
         write(13, "(5ES16.8E2)") (mo_coeff(i), i=e+1, e+nb)
         exit
      end if
   end do


print*, "Calculo realizado con exito"                                           !*** End of option 1


!*******************************Option2: Coordinates Scan**********************************


else if(option1==2) then                                                        !*** If we select option 2
   energy_scan(1:steps)=0
   t=0
   electrons=electrons/2

!*******************************Distance anlysis:

   if (coordenada=="Enlace") then   
      distance(atom_ref)=distance(atom_ref)+real(step_min)*delta
      inicio=distance(atom_ref)
                                                                                
                                                                                !*** lapack subroutines
      allocate(m1(nf, nf), m2(nf, nf), energy(nf), WORK(LWORK))

      do i=1, steps                                                             !*** loop that change the coordinate analysed
         distance(atom_ref)=distance(atom_ref)+real(t)*delta
         t=1
         call z_translate(distance, ang, die, union, r, nro)
         r(1:3, 1:nro)=r(1:3, 1:nro)/a

         call S_H_matrix (S, ham, simb, g, alfa, d, r, iden, norm, nf)

            m1(1:nf, 1:nf)=ham(1:nf, 1:nf)
            m2(1:nf, 1:nf)=S(1:nf, 1:nf)
         LWORK=-1
         call dsygv (1, 'V', 'L', nf, m1, nf, m2, nf, energy, WORK, LWORK, INFO)
         LWORK=WORK(1)
         deallocate (WORK)
         allocate (WORK(LWORK))
         call dsygv (1, 'V', 'L', nf, m1, nf, m2, nf, energy, WORK, LWORK, INFO)

         if (INFO/=0) then
            print*, "Error during calculation of eigenvectors"
            stop
         end if
         
         do j=1, electrons
            energy_scan(i)=energy_scan(i)+2*energy(j)
         end do
      end do                                                                    !*** end of enlace loop


!*******************************Angle anlysis:    
  
   else if(coordenada=="Angulo") then                                           
      delta=delta*0.01745329252
      ang(atom_ref)=ang(atom_ref)+real(step_min)*delta
      inicio=ang(atom_ref)/0.01745329252
                                                                                 !*** lapack subroutines
      allocate(m1(nf, nf), m2(nf, nf), energy(nf), WORK(LWORK))

                                                                                 !*** loop that change the coordinate analysed
      do i=1, steps
         ang(atom_ref)=ang(atom_ref)+real(t)*delta
         t=1
         call z_translate(distance, ang, die, union, r, nro)
         r(1:3, 1:nro)=r(1:3, 1:nro)/a
         call S_H_matrix (S, ham, simb, g, alfa, d, r, iden, norm, nf)
            m1(1:nf, 1:nf)=ham(1:nf, 1:nf)
            m2(1:nf, 1:nf)=S(1:nf, 1:nf)
         LWORK=-1
         call dsygv (1, 'V', 'L', nf, m1, nf, m2, nf, energy, WORK, LWORK, INFO)
         LWORK=WORK(1)
         deallocate (WORK)
         allocate (WORK(LWORK))
         call dsygv (1, 'V', 'L', nf, m1, nf, m2, nf, energy, WORK, LWORK, INFO)

         if (INFO/=0) then
            print*, "Error during calculation of eigenvectors"
            stop
         end if

         do j=1, electrons
            energy_scan(i)=energy_scan(i)+2*energy(j)
         end do
      end do                                                                     !*** End of angle loop

   
!*******************************Dihedral anlysis:

   else if(coordenada=="Diedro") then                                            
      delta=delta*0.01745329252
      do i=1, nro
         if (union(i, 2)==union(atom_ref, 2)) then
            die(i)=die(i)+real(step_min)*delta
         end if
      end do
      inicio=die(atom_ref)/0.01745329252
                                                                                  !*** lapack subroutines
      allocate(m1(nf, nf), m2(nf, nf), energy(nf), WORK(LWORK))
                                                                                  
      do i=1, steps                                                               !*** loop that change the coordinate analysed
         do j=1, nro
            if (union(j, 2)==union(atom_ref, 2)) then
               die(j)=die(j)+real(t)*delta
               t=1
            end if
         end do

         call z_translate(distance, ang, die, union, r, nro)
         r(1:3, 1:nro)=r(1:3, 1:nro)/a
         call S_H_matrix (S, ham, simb, g, alfa, d, r, iden, norm, nf)
            m1(1:nf, 1:nf)=ham(1:nf, 1:nf)
            m2(1:nf, 1:nf)=S(1:nf, 1:nf)
         LWORK=-1
         call dsygv (1, 'V', 'L', nf, m1, nf, m2, nf, energy, WORK, LWORK, INFO)
         LWORK=WORK(1)
         deallocate (WORK)
         allocate (WORK(LWORK))
         call dsygv (1, 'V', 'L', nf, m1, nf, m2, nf, energy, WORK, LWORK, INFO)

         if (INFO/=0) then
            print*, "Error during calculation of eigenvectors"
            stop
         end if
 
         do j=1, electrons
            energy_scan(i)=energy_scan(i)+2*energy(j)
         end do
      end do                                                                      !*** End of angle loop
   end if


!**********************Creating output file:

   doc=indoc //"-scan"//extout
   doc=adjustl(doc)

   open(unit=11, file=doc, status="new", action="write", iostat=ios)
      if (ios/=0) then
         print*, "Error to create:"
         print*, doc
         stop
      end if


   write(11, "(A)") "Gracias por usar TUCKEL"
   write(11, "(A)") "----------------------------------------------------------"
   write(11, "(A)") "TUCKEL  Copyright (C) 2015  Bustamante Carlos Mauricio"
   write(11, "(A)") "----------------------------------------------------------"
   write(11, "(A)")
   write(11, "(A)") "**********************************************************"
   write(11, "(A)") filename
   write(11, "(A)") "**********************************************************"
   write(11, "(A)") "Analisis de la Energía en funcion de la coordenada elegida"
   write(11, "(A)") "**********************************************************"
   if (coordenada=="Enlace") then
      write(11, "(A)") "Longitud de enlace del átomo "//simb(atom_ref)//"("//trim(str(atom_ref))//")"
      write(11, "(A)") "*******************************************************"
      write(11, "(A)") "Coordenada  E(eV)"
   else if(coordenada=="Angulo") then
      delta=delta/0.01745329252
      write(11, "(A)") "Angulo de enlace del átomo "//simb(atom_ref)//"("//trim(str(atom_ref))//")"
      write(11, "(A)") "*******************************************************"
      write(11, "(A)") "Coordinate  E(eV)"
   else if(coordenada=="Diedro") then
      delta=delta/0.01745329252
      write(11, "(A)") "Angulo diedro del átomo "//simb(atom_ref)//"("//trim(str(atom_ref))//")"
      write(11, "(A)") "*******************************************************"
      write(11, "(A)") "Coordinate  E(eV)"        
   end if

   do i=1, steps
      write(11, "(F10.3, F10.3)") inicio+real(i-1)*delta, energy_scan(i)
   end do

   print*, "Calculo realizado con exito"

end if !***********************************************End of options

!************************Closing units and deallocating variables************************************

close(unit=10)
close(unit=11)
close(unit=13)

deallocate(ai, di, alfa, d, r, gi, g, iden, S, ham, WORK, energy, m1, m2, nval, simb)

if (option1==2) then
   deallocate(energy_scan, distance, ang, die, union)
else if (option2==2) then
   deallocate(distance, ang, die, union)
else if (option1==1) then
   deallocate(P, orden, z_atom, r_cc, shell_type,shell_iden, prim_exp, contr_coeff, mo_coeff)
end if


end program TUCKEL                         !*** End of program



!************** Functions that calculate <gi|gj> integrals:

function s2(x1, y1, z1, x2, y2, z2, a1, a2, d1, d2)
   implicit none
   real*8, intent(in):: x1, y1, z1, x2, y2, z2, a1, a2, d1, d2
   real*8:: s2, ex
   
   ex=exp(-((a1*a2)/(a1+a2))*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))
   s2=(d1*d2*(2**(1.5)*(a1*a2)**(0.75)*(a1+a2)**(-1.5)*ex))

end function s2

function sp (x1, y1, z1, x2, y2, z2, a1, a2, d1, d2)
   implicit none
   real*8, intent(in):: x1, y1, z1, x2, y2, z2, a1, a2, d1, d2
   real*8:: sp, ex
   
   ex=exp(-((a1*a2)/(a1+a2))*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))
   sp= ((d1*d2*(4*2**0.5*a1**1.75*a2**1.25*ex*(x1-x2)))*((a1+a2)**(-2.5)))

end function sp

function p2 (x1, y1, z1, x2, y2, z2, a1, a2, d1, d2)
   implicit none
   real*8, intent(in):: x1, y1, z1, x2, y2, z2, a1, a2, d1, d2
   real*8:: p2, ex

   ex=exp(-((a1*a2)/(a1+a2))*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))   
   p2= (d1*d2*((4*2**0.5*(a1*a2)**1.25*ex*(a2+a1*(1-2*a2*(x1-x2)**2))*((a1+a2)**(-3.5)))))

end function p2

function pp (x1, y1, z1, x2, y2, z2, a1, a2, d1, d2)
   implicit none
   real*8, intent(in):: x1, y1, z1, x2, y2, z2, a1, a2, d1, d2
   real*8:: pp, ex

   ex=exp(-((a1*a2)/(a1+a2))*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))
   pp= (d1*d2*((-8*(2**0.5)*((a1*a2)**2.25)*ex*(x1-x2)*(y1-y2))*((a1+a2)**(-3.5))))


end function pp

!**********Subroutine that order the parameters of Z-matrix:

subroutine z_elements(distance, simb, ang, die, union, nro)
      implicit none
      real*8, allocatable, intent(inout):: distance(:), ang(:), die(:)
      integer, allocatable, intent(inout):: union(:, :)
      character(len=4), allocatable, intent(inout):: simb(:)   
      integer, intent(in):: nro
      integer :: i

   union(1, 1)=1
   read(10, *), simb(1)
   union(2, 1)=2
   read(10, *), simb(2), union(2, 2), distance(2)
   union(3, 1)=3
   read(10, *), simb(3), union(3, 2), distance(3), union(3, 3), ang(3)

   do i=4, nro
      union(i, 1)=i
      read(10, *), simb(i), union(i, 2), distance(i), union(i, 3), ang(i), union(i, 4), die(i)
      die(i)=-die(i)
   end do

   ang(1:nro)=ang(1:nro)*0.01745329252
   die(1:nro)=die(1:nro)*0.01745329252

end subroutine z_elements


!***** Subroutine that translate the Z-coordinates into cartesian coordinates:
subroutine z_translate(distance, ang, die, union, r, nro)
   
   implicit none
   
   real*8, allocatable, intent(in):: distance(:), ang(:), die(:)
   real*8, allocatable, intent(inout):: r(:, :)
   real*8:: v1(3), v2(3), v3(3), r1(3), r2(3), r3(3), x, y, z, nule=0
   integer, allocatable, intent(in):: union(:, :)
   integer, intent(in):: nro
   integer:: i
   

   r(1:3, 1:nro)=0
   r(1:3, 2)=(/nule, nule, distance(2)/)

      if (union(3, 2)==1)then
         r(1:3, 3)=(/(distance(3)*sin(ang(3))), nule, (distance(3)*cos(ang(3)))/)
      else
         r(1:3, 3)=(/(distance(3)*sin(ang(3))), nule, -(distance(3)*cos(ang(3)))+r(3, 2)/)
      end if

      do i=4, nro
         r(1:3, i)=(/distance(i)*sin(ang(i))*cos(die(i)), distance(i)*sin(ang(i))*sin(die(i)), distance(i)*cos(ang(i))/)
      end do

      do i=4, nro
      !   if(union(i, 2)/=1)then
            r1(1:3)=r(1:3, union(i, 3))
            r2(1:3)=r(1:3, union(i, 2))
            r3(1:3)=r(1:3, union(i, 4))
            call unitvector(r1, r2, r3, v1, v2, v3)
            x=r(1, i)*(v1(1))+r(2, i)*(v2(1))+r(3, i)*(v3(1))
            y=r(1, i)*(v1(2))+r(2, i)*(v2(2))+r(3, i)*(v3(2))
            z=r(1, i)*(v1(3))+r(2, i)*(v2(3))+r(3, i)*(v3(3))
            r(1:3, i)=(/(x+r(1, union(i, 2))), (y+r(2, union(i, 2))), (z+r(3, union(i, 2)))/)
       !  end if
      end do

end subroutine z_translate


!*** Subroutine that create S array and finish H array:

subroutine S_H_matrix (S, ham, simb, g, alfa, d, r, iden, norm, nf)

   implicit none
   
   real*8, allocatable, intent(inout):: S(:, :), norm(:), ham(:, :)
   real*8, allocatable, intent(in):: d(:, :), alfa(:, :), r(:, :)
   character (len=4), allocatable, intent(in):: g(:, :), simb(:)
   integer, allocatable, intent(in):: iden(:)
   integer, intent(in):: nf
   integer:: i, j, k, l
   real*8:: s2, sp, p2, pp
   real*8:: x1, y1, z1, x2, y2, z2, a1, a2, d1, d2, total
   

   S(1:nf, 1:nf)=0
   norm(1:nf)=0

     do i=1, nf
           x1=r(1, iden(i))
           y1=r(2, iden(i))
           z1=r(3, iden(i))

           do k=1, 6
               do l=1, 6
                  a1=alfa(k, i)
                  a2=alfa(l, i)
                  d1=d(k, i)
                  d2=d(l, i)
  
                  if (g(k, i)=="ss") then
                      total= s2(x1, y1, z1, x1, y1, z1, a1, a2, d1, d2)
                  else if (g(k, i)=="px") then
                      total= p2(x1, y1, z1, x1, y1, z1, a1, a2, d1, d2)
                  else if (g(k, i)=="py") then
                      total= p2(y1, x1, z1, y1, x1, z1, a1, a2, d1, d2)
                  else if (g(k, i)=="pz") then
                      total= p2(z1, y1, x1, z1, y1, x1, a1, a2, d1, d2)
                  end if
                      S(i, i)=S(i, i)+ total
               end do
           end do
              norm(i)=(S(i, i)**(-0.5))
     end do

   S(1:nf, 1:nf)=0
      do i=1, nf
         do j=i, nf
            x1=r(1, iden(i))
            y1=r(2, iden(i))
            z1=r(3, iden(i))
            x2=r(1, iden(j))
            y2=r(2, iden(j))
            z2=r(3, iden(j))

            do k=1, 6
                do l=1, 6
                  a1=alfa(k, i)
                  a2=alfa(l, j)
                  d1=d(k, i)
                  d2=d(l, j)

               if (g(k, i)=="ss") then

                  if (g(l, j)=="ss") then
                     total= s2(x1, y1, z1, x2, y2, z2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  else if (g(l, j)=="px") then
                     total= sp(x1, y1, z1, x2, y2, z2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  else if (g(l, j)=="py") then
                     total= sp(y1, x1, z1, y2, x2, z2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  else if (g(l, j)=="pz") then
                     total= sp(z1, y1, x1, z2, y2, x2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  end if
               end if
               if (g(k, i)=="px") then
                  if (g(l, j)=="ss") then
                     total= sp(x2, y2, z2, x1, y1, z1, a2, a1, d2, d1)*(norm(i)*norm(j))
                  else if (g(l, j)=="px") then
                     total= p2(x1, y1, z1, x2, y2, z2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  else if (g(l, j)=="py") then
                     total= pp(x1, y1, z1, x2, y2, z2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  else if (g(l, j)=="pz") then
                     total= pp(x1, z1, y1, x2, z2, y2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  end if
               end if
               if (g(k, i)=="py") then
                  if (g(l, j)=="ss") then
                     total= sp(y2, x2, z2, y1, x1, z1, a2, a1, d2, d1)*(norm(i)*norm(j))
                  else if (g(l, j)=="px") then
                     total= pp(y1, x1, z1, y2, x2, z2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  else if (g(l, j)=="py") then
                     total= p2(y1, x1, z1, y2, x2, z2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  else if (g(l, j)=="pz") then
                     total= pp(y1, z1, x1, y2, z2, x2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  end if
               end if
               if (g(k, i)=="pz") then
                  if (g(l, j)=="ss") then
                     total= sp(z2, y2, x2, z1, y1, x1, a2, a1, d2, d1)*(norm(i)*norm(j))
                  else if (g(l, j)=="px") then
                     total= pp(z1, x1, y1, z2, x2, y2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  else if (g(l, j)=="py") then
                     total= pp(z1, y1, x1, z2, y2, x2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  else if (g(l, j)=="pz") then
                     total= p2(z1, y1, x1, z2, y2, x2, a1, a2, d1, d2)*(norm(i)*norm(j))
                  end if
               end if

                  S(i, j)=S(i, j)+ total

               end do
            end do
         end do
      end do

   do i=1, nf-1
      do j=i+1, nf
         S(j, i)= S(i, j)
      end do
   end do

   do i=1, nf-1
      do j=i+1, nf
         if (simb(i)=="N" .or. simb(j)=="N") then
            ham(i, j)=-35.0*S(i, j)
            ham(j, i)=ham(i, j)
         else
            ham(i, j)=0.5*1.75*(ham(i, i)+ham(j, j))*S(i, j)
            ham(j, i)=ham(i, j)
         end if
      end do
   end do

end subroutine S_H_matrix


!*** Subroutine that assign to each atom the corresponding basis parameters:

subroutine element (atom, gi, ai, di, orbitals)
   !STO-6G
   implicit none

   character (len=1), intent(in) :: atom
   character (len=4) :: C="C", H="H", B="B", N="N", O="O", F="F"
   character(len=4), allocatable, intent(out):: gi(:, :)
   real*8, allocatable, intent(out):: ai(:, :), di(:, :)
   integer, intent(out):: orbitals

   if (atom == C)then
      allocate (ai(6, 4), di(6, 4), gi(6, 4))
      orbitals= 4

      gi(1, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(2, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(3, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(4, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(5, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(6, 1:4)=(/"ss", "px", "py", "pz"/)

      !sofi ama a charly
      ai(1, 1:4)=(/30.49723950, 30.49723950, 30.49723950, 30.49723950/)
      ai(2, 1:4)=(/6.036199601, 6.036199601, 6.036199601, 6.036199601/)
      ai(3, 1:4)=(/1.876046337, 1.876046337, 1.876046337, 1.876046337/)
      ai(4, 1:4)=(/0.7217826470, 0.7217826470, 0.7217826470, 0.7217826470/)
      ai(5, 1:4)=(/0.3134706954, 0.3134706954, 0.3134706954, 0.3134706954/)
      ai(6, 1:4)=(/0.1436865550, 0.1436865550, 0.1436865550, 0.1436865550/)

      di(1, 1:4)=(/-0.01325278809, 0.0037596966, 0.0037596966, 0.0037596966/)
      di(2, 1:4)=(/-0.04699171014, 0.0376793698, 0.0376793698, 0.0376793698/)
      di(3, 1:4)=(/-0.03378537151, 0.1738967435, 0.1738967435, 0.1738967435/)
      di(4, 1:4)=(/0.25024178610, 0.4180364347, 0.4180364347, 0.4180364347/)
      di(5, 1:4)=(/0.59511725260, 0.4258595477, 0.4258595477, 0.4258595477/)
      di(6, 1:4)=(/0.24070617630, 0.1017082955, 0.1017082955, 0.1017082955/)

   else if (atom == H) then
      orbitals= 1
      allocate(ai(6, 1), di(6, 1), gi(6, 1))

      gi(1, 1)="ss"
      gi(2, 1)="ss"
      gi(3, 1)="ss"
      gi(4, 1)="ss"
      gi(5, 1)="ss"
      gi(6, 1)="ss"

      ai(1, 1)=35.52322122
      ai(2, 1)=6.513143725 
      ai(3, 1)=1.822142904 
      ai(4, 1)=0.625955266
      ai(5, 1)=0.243076747
      ai(6, 1)=0.100112428

      di(1, 1)=0.00916359628
      di(2, 1)=0.04936149294 
      di(3, 1)=0.16853830490
      di(4, 1)=0.37056279970
      di(5, 1)=0.41649152980
      di(6, 1)=0.13033408410 

   else if (atom == B)then
      allocate (ai(6, 4), di(6, 4), gi(6, 4))
      orbitals= 4

      gi(1, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(2, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(3, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(4, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(5, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(6, 1:4)=(/"ss", "px", "py", "pz"/)

      ai(1, 1:4)=(/23.194575000, 23.194575000, 23.194575000, 23.194575000/)
      ai(2, 1:4)=(/4.5908100000, 4.5908100000, 4.5908100000, 4.5908100000/)
      ai(3, 1:4)=(/1.4268195000, 1.4268195000, 1.4268195000, 1.4268195000/)
      ai(4, 1:4)=(/0.5489482500, 0.5489482500, 0.5489482500, 0.5489482500/)
      ai(5, 1:4)=(/0.2384100000, 0.2384100000, 0.2384100000, 0.2384100000/)
      ai(6, 1:4)=(/0.1092802500, 0.1092802500, 0.1092802500, 0.1092802500/)

      di(1, 1:4)=(/-0.01325278809, 0.00375969662, 0.00375969662, 0.00375969662/)
      di(2, 1:4)=(/-0.04699171014, 0.03767936984, 0.03767936984, 0.03767936984/)
      di(3, 1:4)=(/-0.03378537151, 0.17389674350, 0.17389674350, 0.17389674350/)
      di(4, 1:4)=(/0.25024178610, 0.41803643470, 0.41803643470, 0.41803643470/)
      di(5, 1:4)=(/0.59511725260, 0.42585954770, 0.42585954770, 0.42585954770/)
      di(6, 1:4)=(/0.24070617630, 0.10170829550, 0.10170829550, 0.10170829550/)

   else if (atom == N)then
      allocate (ai(6, 4), di(6, 4), gi(6, 4))
      orbitals= 4

      gi(1, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(2, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(3, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(4, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(5, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(6, 1:4)=(/"ss", "px", "py", "pz"/)

      ai(1, 1:4)=(/39.19880787, 39.19880787, 39.19880787, 39.19880787/)
      ai(2, 1:4)=(/7.758467071, 7.758467071, 7.758467071, 7.758467071/)
      ai(3, 1:4)=(/2.411325783, 2.411325783, 2.411325783, 2.411325783/)
      ai(4, 1:4)=(/0.9277239437, 0.9277239437, 0.9277239437, 0.9277239437/)
      ai(5, 1:4)=(/0.4029111410, 0.4029111410, 0.4029111410, 0.4029111410/)
      ai(6, 1:4)=(/0.1846836552, 0.1846836552, 0.1846836552, 0.1846836552/)

      di(1, 1:4)=(/-0.01325278809, 0.0037596966, 0.0037596966, 0.0037596966/)
      di(2, 1:4)=(/-0.04699171014, 0.0376793698, 0.0376793698, 0.0376793698/)
      di(3, 1:4)=(/-0.03378537151, 0.1738967435, 0.1738967435, 0.1738967435/)
      di(4, 1:4)=(/0.25024178610, 0.4180364347, 0.4180364347, 0.4180364347/)
      di(5, 1:4)=(/0.59511725260, 0.4258595477, 0.4258595477, 0.4258595477/)
      di(6, 1:4)=(/0.24070617630, 0.1017082955, 0.1017082955, 0.1017082955/)

   else if (atom == O)then
      allocate (ai(6, 4), di(6, 4), gi(6, 4))
      orbitals= 4

      gi(1, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(2, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(3, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(4, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(5, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(6, 1:4)=(/"ss", "px", "py", "pz"/)

      ai(1, 1:4)=(/52.18776196, 52.18776196, 52.18776196, 52.18776196/)
      ai(2, 1:4)=(/10.32932006, 10.32932006, 10.32932006, 10.32932006/)
      ai(3, 1:4)=(/3.210344977, 3.210344977, 3.210344977, 3.210344977/)
      ai(4, 1:4)=(/1.235135428, 1.235135428, 1.235135428, 1.235135428/)
      ai(5, 1:4)=(/0.536420158, 0.536420158, 0.536420158, 0.536420158/)
      ai(6, 1:4)=(/0.245880606, 0.245880606, 0.245880606, 0.245880606/)

      di(1, 1:4)=(/-0.01325278809, 0.0037596966, 0.0037596966, 0.0037596966/)
      di(2, 1:4)=(/-0.04699171014, 0.0376793698, 0.0376793698, 0.0376793698/)
      di(3, 1:4)=(/-0.03378537151, 0.1738967435, 0.1738967435, 0.1738967435/)
      di(4, 1:4)=(/0.25024178610, 0.4180364347, 0.4180364347, 0.4180364347/)
      di(5, 1:4)=(/0.59511725260, 0.4258595477, 0.4258595477, 0.4258595477/)
      di(6, 1:4)=(/0.24070617630, 0.1017082955, 0.1017082955, 0.1017082955/)

   else if (atom == F)then
      allocate (ai(6, 4), di(6, 4), gi(6, 4))
      orbitals= 4

      gi(1, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(2, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(3, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(4, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(5, 1:4)=(/"ss", "px", "py", "pz"/)
      gi(6, 1:4)=(/"ss", "px", "py", "pz"/)

      ai(1, 1:4)=(/67.03228091, 67.03228091, 67.03228091, 67.03228091/)
      ai(2, 1:4)=(/13.26743777, 13.26743777, 13.26743777, 13.26743777/)
      ai(3, 1:4)=(/4.123509771, 4.123509771, 4.123509771, 4.123509771/)
      ai(4, 1:4)=(/1.586462839, 1.586462839, 1.586462839, 1.586462839/)
      ai(5, 1:4)=(/0.689001892, 0.689001892, 0.689001892, 0.689001892/)
      ai(6, 1:4)=(/0.315819978, 0.315819978, 0.315819978, 0.315819978/)

      di(1, 1:4)=(/-0.01325278809, 0.0037596966, 0.0037596966, 0.0037596966/)
      di(2, 1:4)=(/-0.04699171014, 0.0376793698, 0.0376793698, 0.0376793698/)
      di(3, 1:4)=(/-0.03378537151, 0.1738967435, 0.1738967435, 0.1738967435/)
      di(4, 1:4)=(/0.25024178610, 0.4180364347, 0.4180364347, 0.4180364347/)
      di(5, 1:4)=(/0.59511725260, 0.4258595477, 0.4258595477, 0.4258595477/)
      di(6, 1:4)=(/0.24070617630, 0.1017082955, 0.1017082955, 0.1017082955/)

   end if
return
end subroutine element


!*** function that translate integers into characters:

function str(t)
    integer, intent(in) :: t
    character(len=20):: str
    write (str, *) t
    str = adjustl(str)
end function str


!*** subroutine that creates unitary vectors:

subroutine unitvector (r1, r2, r3, v1, v2, v3)

   implicit none
   real*8, intent(out):: v1(3), v2(3), v3(3)
   real*8:: r12(3), r13(3), pro
   real*8, intent(in):: r1(3), r2(3), r3(3)

   r12(1:3)=r1(1:3)-r2(1:3)
   r13(1:3)=r3(1:3)-r1(1:3)

   call product_x(r12, r13, v2)

   v3(1:3)=r12(1:3)/((r12(1)**2+r12(2)**2+r12(3)**2)**0.5)

   call product_x(v2, v3, v1)

   pro=r13(1)*v1(1)+r13(2)*v1(2)+r13(3)*v1(3)
   if (pro<0)then
      v1=-v1
      v2=-v2
   end if

end subroutine unitvector

!*** Subroutine that solve vectorial products:
subroutine product_x (a, b, n)
        implicit none
        real*8, dimension(3), intent(in):: a, b
        real*8, dimension(3), intent(out):: n
        real*8::norm

        n(1)=a(2)*b(3)-a(3)*b(2)
        n(2)=a(3)*b(1)-a(1)*b(3)
        n(3)=a(1)*b(2)-a(2)*b(1)

        norm=(n(1)**2+n(2)**2+n(3)**2)**0.5
        n(1:3)=n(1:3)/norm

end subroutine product_x
