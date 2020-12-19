#	Hasan Zerze November 17, 2015
#	This code locates atoms on the surface of a colloid
#	Given atom separation (d) and the radius of the colloid (Rc)
#		By only changing the value of "Rc" and keeping "d" constant, one can get constant surface density for different colloid sizes
#	Use: ./colloidsurface < input.dat
#	Input file example:
#	1.77 5.0
#	for a case in which d=1.77 and Rc=5.0
#	Output files are colloid.pdb and colloid.xyz
#	unitsphere.pdb and unitsphere.xyz files represent unit radius form of the same colloid
#	So, all coordinates in the unitsphere outputs were divided by the colloid radius to normalize
#	
	subroutine colloid(d, Rc, totalnum)
	double precision d, Rc, r
        double precision theta, dtheta, phi, dphi
	double precision x, y, z
        double precision PI
	character *50 frmt1, frmt2, frmt3
	integer i,j, Ni, Nj, iatom,natom
	integer totalnum
Cf2py	intent(out) totalnum

	PI = 4.0d0 * atan(1.0d0)
	open(10, file='colloid.pdb')
        open(11, file='colloid.xyz')
	open(12, file='unitsphere.pdb')
        open(13, file='unitsphere.xyz')
	frmt1="(A4,I7,A3,A8,A4,F12.3, 2F8.3, 2F6.2,A12)"

C	read(*,*) d, Rc

CCC	Calculating number of latitudes (Ni):

	dtheta = d / Rc
	
	Ni = nint( PI / dtheta ) + 1

	dtheta = PI / ( Ni - 1 ) 	! fine-tuned to make all particles equally spaced

CCC	loop over all the latitudes:
	iatom = 1
        write(10,*) 'CRYST1    0.000    0.000    0.000  
     &          90.00  90.00  90.00 P 1           1'
        write(12,*) 'CRYST1    0.000    0.000    0.000
     &          90.00  90.00  90.00 P 1           1'

	do i = 0, Ni-1

		theta = (PI / 2.0) - i * dtheta
		r = Rc * cos(theta)

CCC	Determine the number of atoms at the present latitude:

		if(i.ne.0.and.i.ne.(Ni-1)) then		! means except for the pole positions
			dphi = d / r
			Nj = nint(2.0d0 * PI / dphi) + 1
			dphi = 2.0d0 * PI / (Nj - 1)	! fine-tuned to make all particles equally spaced
		else
			Nj = 1
		endif

CCC	Now loop over all the atoms for the present latitude
		do j = 0, Nj-1
			phi = j * dphi
			x = r * cos(phi)
			y = Rc * sin(theta)
			z = r * sin(phi)
             write(10,frmt1) 'ATOM', iatom, 'H','X','1'       
     &		, x, y, z,0.00,0.00,'H'
             write(12,frmt1) 'ATOM', iatom, 'H','X','1'
     &          , x/Rc, y/Rc, z/Rc,0.00,0.00,'H'

			iatom = iatom + 1
		enddo

	enddo

        write(10,*) 'END'
        write(12,*) 'END'

	natom=iatom-1
	totalnum=natom

CC	Same algorithm for xyz file
	iatom=0
        write(11,*) natom
        write(11,*) 'Atoms'
        write(13,*) natom
        write(13,*) 'Atoms'

        do i = 0, Ni-1

                theta = (PI / 2.0) - i * dtheta
                r = Rc * cos(theta)

CCC     Determine the number of atoms at the present latitude:

                if(i.ne.0.and.i.ne.(Ni-1)) then         ! means except for the pole positions
                        dphi = d / r
                        Nj = nint(2.0d0 * PI / dphi) + 1
                        dphi = 2.0d0 * PI / (Nj - 1)    ! fine-tuned to make all particles equally spaced
                else
                        Nj = 1
                endif

CCC     Now loop over all the atoms for the present latitude
                do j = 0, Nj-1
                        phi = j * dphi
                        x = r * cos(phi)
                        y = Rc * sin(theta)
                        z = r * sin(phi)
                        write(11,*) '1', x, y, z
                        write(13,*) '1', x/Rc, y/Rc, z/Rc

                        iatom = iatom + 1
                enddo

        enddo
	
	close (10)
        close (11)
        close (12)
        close (13)	

	end
