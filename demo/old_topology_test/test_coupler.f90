module test_inputs
implicit none

	integer :: npx_md 
	integer :: npy_md 
	integer :: npz_md 
	integer :: npx_cfd 
	integer :: npy_cfd 
	integer :: npz_cfd 
	integer :: ncx 
	integer :: ncy 
	integer :: ncz 
	real(kind(0.d0)) :: xL_md 
	real(kind(0.d0)) :: yL_md 
	real(kind(0.d0)) :: zL_md
	real(kind(0.d0)) :: xL_cfd 
	real(kind(0.d0)) :: yL_cfd 
	real(kind(0.d0)) :: zL_cfd

end module test_inputs

!---- TEST PROGRAM -----------------!
program test_coupler
	implicit none

	call initialise
	call test_setup
	call test_gather_scatter
	call test_send_recv_MD2CFD
	call test_send_recv_CFD2MD
	call finalise

end program test_coupler
!-----------------------------------!

subroutine test_setup
	use mpi
	use test_inputs
	implicit none

	integer :: ierr
	integer :: myid_test, rank_test

	open(unit=1,file='TOPOL.in',form='formatted')
		read(1,*) npx_md 
		read(1,*) npy_md 
		read(1,*) npz_md 
		read(1,*) npx_cfd 
		read(1,*) npy_cfd 
		read(1,*) npz_cfd 
		read(1,*) ncx 
		read(1,*) ncy 
		read(1,*) ncz 
		read(1,*) xL_md 
		read(1,*) yL_md 
		read(1,*) zL_md 
		read(1,*) xL_cfd 
		read(1,*) yL_cfd 
		read(1,*) zL_cfd 
	close(1)

	call MPI_COMM_RANK(MPI_COMM_WORLD,myid_test,ierr)
	rank_test = myid_test + 1

	if (rank_test .le. npx_cfd*npy_cfd*npz_cfd) then
		call test_cfd_setup	
	else
		call test_md_setup
	end if

contains

	subroutine test_cfd_setup
	use CPL, only: CPL_create_comm, cfd_realm, coupler_cfd_init
	implicit none

		! Local vars
		integer :: nproc_cfd
		integer :: realm_comm 
		integer :: cart_comm,trank,tid,tcoord(3),ngx,ngy,ngz
		integer :: dims(3),periods(3)
		integer :: i,j,k
		real(kind(0.d0)) :: dx, dy, dz

		! CFD init inputs
		integer :: nsteps,icomm_grid 
		integer,dimension(3) :: ijkcmin,ijkcmax,npxyz_cfd,ncxyz
		integer,dimension(:), allocatable :: iTmin,iTmax,jTmin,jTmax,kTmin,kTmax
		integer,dimension(:,:), allocatable :: icoord
		real(kind(0.d0)) :: dt,density
		real(kind(0.d0)),dimension(3) :: xyzL
		real(kind(0.d0)),dimension(:), allocatable :: zgrid
		real(kind(0.d0)),dimension(:,:), allocatable :: xgrid,ygrid

		call CPL_create_comm(cfd_realm,realm_comm,ierr)
	
		nproc_cfd = npx_cfd * npy_cfd * npz_cfd

		dims = (/npx_cfd, npy_cfd, npz_cfd/)
		periods = (/1, 0, 1/)
		call MPI_CART_CREATE(realm_comm,3,dims,periods,.true.,cart_comm,ierr)

		! Setup coupler_cfd_init inputs
		nsteps = 1
		dt = 0.1d0
		icomm_grid = cart_comm 
		allocate(icoord(3,nproc_cfd))	
		do trank = 1,nproc_cfd
			tid = trank - 1	
			call MPI_CART_COORDS(cart_comm,tid,3,tcoord,ierr)
			tcoord = tcoord + 1
			icoord(:,trank) = tcoord(:)
		end do 
		npxyz_cfd = (/npx_cfd,npy_cfd,npz_cfd/)
		xyzL = (/xL_cfd,yL_cfd,zL_cfd/)
		ncxyz = (/ncx,ncy,ncz/)
		density = 0.8 
		ijkcmax = (/ncx,ncy,ncz/) 
		ijkcmin = (/1,1,1/) 
		allocate(iTmin(npx_cfd))
		allocate(iTmax(npx_cfd))
		allocate(jTmin(npy_cfd))
		allocate(jTmax(npy_cfd))
		allocate(kTmin(npz_cfd))
		allocate(kTmax(npz_cfd))
		do i = 1, npx_cfd
			iTmax(i) = i * (ncx / npx_cfd)
			iTmin(i) = iTmax(i) - (ncx / npx_cfd) + 1 
		end do
		do j = 1, npy_cfd
			jTmax(j) = j * (ncy / npy_cfd)
			jTmin(j) = jTmax(j) - (ncy / npy_cfd) + 1 
		end do
		do k = 1, npz_cfd
			kTmax(k) = k * (ncz / npz_cfd)
			kTmin(k) = kTmax(k) - (ncz / npz_cfd) + 1 
		end do
		ngx = ncx + 1
		ngy = ncy + 1
		ngz = ncz + 1
		dx = xL_cfd / real(ncx,kind(0.d0))
		dy = yL_cfd / real(ncy,kind(0.d0))
		dz = zL_cfd / real(ncz,kind(0.d0))
		allocate(xgrid(ngx,ngy))
		allocate(ygrid(ngx,ngy))
		allocate(zgrid(ngz))
		do i = 1,ngx
		do j = 1,ngy
			xgrid(i,j) = ( i - 1 ) * dx 
			ygrid(i,j) = ( j - 1 ) * dy
		end do
		end do
		do k = 1,ngz
			zgrid( k ) = ( k - 1 ) * dz
		end do
	
		call coupler_cfd_init(nsteps,dt,icomm_grid,icoord,npxyz_cfd,xyzL,ncxyz, & 
		                      density,ijkcmax,ijkcmin,iTmin,iTmax,jTmin, & 
		                      jTmax,kTmin,kTmax,xgrid,ygrid,zgrid)

	end subroutine test_cfd_setup

	subroutine test_md_setup
	use CPL, only: CPL_create_comm, md_realm, coupler_md_init
	implicit none

		! Local vars
		integer :: nproc_md
		integer :: realm_comm
		integer :: cart_comm, trank, tid, tcoord(3)
		integer :: dims(3), periods(3)

		! MD init inputs
		integer :: nsteps, icomm_grid,initialstep
		integer,dimension(3) :: npxyz_md	
		integer,dimension(:,:),allocatable :: icoord
		real(kind(0.d0)) :: dt,density
		real(kind=kind(0.d0)),dimension(3) :: globaldomain

		call CPL_create_comm(md_realm,realm_comm,ierr)

		nproc_md  = npx_md  * npy_md  * npz_md	
	
		dims = (/npx_md,npy_md,npz_md/)
		periods = (/1, 0, 1/)
		call MPI_CART_CREATE(realm_comm,3,dims,periods,.true.,cart_comm,ierr)

		! Setup coupler_md_init inputs
		nsteps = 1
		initialstep = 1
		dt = 0.1d0
		icomm_grid = cart_comm
		allocate(icoord(3,nproc_md))	
		do trank = 1,nproc_md
			tid = trank - 1	
			call MPI_CART_COORDS(cart_comm,tid,3,tcoord,ierr)
			tcoord = tcoord + 1
			icoord(:,trank) = tcoord(:)
		end do 
		npxyz_md = (/npx_md,npy_md,npz_md/)
		globaldomain = (/xL_md,yL_md,zL_md/)
		density = 0.8

		call coupler_md_init(nsteps,initialstep,dt,icomm_grid,icoord,npxyz_md, &
		                     globaldomain,density)

!		call write_overlap_comms_md

	end subroutine test_md_setup

end subroutine test_setup

subroutine write_realm_info
	use CPL
	use mpi
	implicit none

	integer :: coord(3)

	if (myid_world.eq.0) then
		write(1000+rank_world,*), '---------- REALM INFORMATION --------------'
		write(1000+rank_world,*), ' wrank  realm  realmrank       cart coords '
		write(1000+rank_world,*), '-------------------------------------------'
	end if

	call MPI_cart_coords(CPL_CART_COMM,myid_cart,3,coord,ierr)
	coord(:) = coord(:) + 1

	write(1000+rank_world,'(3i6,a10,3i5)'),rank_world, realm, rank_realm,'',&
	                                       coord(1), coord(2), coord(3)
	
	if (myid_world.eq.nproc_world) then
		write(1000+rank_world,*), '------------ END REALM INFO ---------------'
		write(1000+rank_world,*), '==========================================='
	end if
	
end subroutine write_realm_info

subroutine write_overlap_comms_md
	use CPL
	use mpi
	implicit none

	integer :: coord(3)

	if (myid_realm.eq.0) then
		write(2000+rank_realm,*),'rank_realm,rank_olap,  mdcoord,'  &
		                        ,'   olap_mask,  CPL_OLAP_COMM'
	end if
	
	call MPI_cart_coords(CPL_CART_COMM,myid_cart,3,coord,ierr)
	coord(:) = coord(:) + 1

	write(2000+rank_realm,'(2i7,a5,3i5,a5,l4,a5,i20)'), &
		rank_realm,rank_olap,'',coord,'',CPL_overlap(), &
		'',CPL_OLAP_COMM

end subroutine write_overlap_comms_md

subroutine test_gather_scatter
	use CPL, only: CPL_overlap,realm,md_realm,cfd_realm,ncx,ncy,ncz, &
	               CPL_CART_COMM,jcmax_olap,rank_cart,ierr,olap_mask, &
	               rank_world,myid_world, &
	               CPL_Cart_coords, CPL_proc_extents, CPL_gather, CPL_scatter
	implicit none

	double precision,dimension(:,:,:,:),allocatable	:: u,stress,gatheru,scatterstress
	integer :: coord(3), extents(6), gatherlims(6), scatterlims(6), npercell
	integer :: pos, ixyz, icell, jcell, kcell
	integer :: ncxl,ncyl,nczl
	integer :: i,j,k
	
 	if (.not.CPL_overlap()) return

	if (realm .eq. md_realm) then	

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)
		npercell = 3
		allocate(u(npercell,extents(1):extents(2), &
		                    extents(3):extents(4), &
		                    extents(5):extents(6)))
		allocate(stress(0,0,0,0))

		! Populate dummy gatherbuf
		pos = 1
		do ixyz = 1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)

			u(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                                      1000*jcell + &
			                                   1000000*kcell
			pos = pos + 1

		end do
		end do
		end do
		end do

	else if (realm .eq. cfd_realm) then	  
		
		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,extents)
		npercell = 9
		allocate(u(0,0,0,0))
		allocate(stress(npercell,extents(1):extents(2), &
		                         extents(3):extents(4), &
		                         extents(5):extents(6)))

		! Populate dummy gatherbuf
		pos = 1
		do ixyz = 1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)

			stress(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                                           1000*jcell + &
			                                        1000000*kcell
			pos = pos + 1

		end do
		end do
		end do
		end do

	endif

	! Allocate test arrays over local domain
	if (realm.eq.cfd_realm) then
		call CPL_cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,extents)
		ncxl = extents(2) - extents(1) + 1
		ncyl = extents(4) - extents(3) + 1
		nczl = extents(6) - extents(5) + 1
		allocate(gatheru(3,ncxl,ncyl,nczl))
		gatheru = 0.d0
	else if (realm.eq.md_realm) then
		call CPL_cart_coords(CPL_CART_COMM,rank_cart,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)
		ncxl = extents(2) - extents(1) + 1
		ncyl = extents(4) - extents(3) + 1
		nczl = extents(6) - extents(5) + 1
		allocate(scatterstress(9,ncxl,ncyl,nczl))
		scatterstress = 0.d0
	end if

	!gatherlims  = (/1,1,1,1,1,1/)
	!scatterlims = (/1,1,1,1,1,1/)
	!================== PERFORM GATHER/SCATTER =============================!	
	gatherlims  = (/1,ncx, jcmax_olap, jcmax_olap , 1, ncz/)
	scatterlims = (/1,ncx, 1, 1, 1,ncz/)
	if (olap_mask(rank_world).eqv..true.) call CPL_gather(u,3,gatherlims,gatheru)
	if (olap_mask(rank_world).eqv..true.) call CPL_scatter(stress,9,scatterlims, &
	                                                      scatterstress)

	! Print results to file
	if (realm.eq.cfd_realm) then

		do ixyz  = 1,size(gatheru,1)
		do icell = 1,size(gatheru,2)
		do jcell = 1,size(gatheru,3)
		do kcell = 1,size(gatheru,4)

			i = icell + extents(1) - 1
			j = jcell + extents(3) - 1
			k = kcell + extents(5) - 1

			if (gatheru(ixyz,icell,jcell,kcell).lt.0.0001) then
				!write(8000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				!	  'gatheru(',0,',',0,',',0,',',0,') =', 0.d0
			else
				write(8000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					  'gatheru(',ixyz,',',i,',',j,',',k,') =', &
					   gatheru(ixyz,icell,jcell,kcell)
			end if

		end do	
		end do	
		end do
		end do

	else if (realm.eq.md_realm) then

		do ixyz  = 1,size(scatterstress,1)
		do icell = 1,size(scatterstress,2)
		do jcell = 1,size(scatterstress,3)
		do kcell = 1,size(scatterstress,4)

			i = icell + extents(1) - 1
			j = jcell + extents(3) - 1
			k = kcell + extents(5) - 1

			if (scatterstress(ixyz,icell,jcell,kcell).lt.0.0001) then
				!write(7000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				!	  'scatterstress(',0,',',0,',',0,',',0,') =', 0.d0
			else
				write(7000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					  'scatterstress(',ixyz,',',i,',',j,',',k,') =', &
					   scatterstress(ixyz,icell,jcell,kcell)
			end if

		end do	
		end do	
		end do
		end do
	
	end if
	
end subroutine test_gather_scatter

subroutine test_packing
	use coupler_module
	use coupler
	implicit none


	integer 										:: ncxl, ncyl, nczl
	integer 										:: coord(3), extents(6)
	integer											:: ixyz, icell, jcell, kcell
	double precision,dimension(:),allocatable		:: outbuf
	double precision,dimension(:,:,:,:),allocatable	:: packbuf,testbuf

	! Test Packing in the overlap region only
	if (olap_mask(rank_world).neqv..true.) return
				   
	if (realm .eq. md_realm) then		   

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)

		ncxl = extents(2)-extents(1)+1
		ncyl = extents(4)-extents(3)+1
		nczl = extents(6)-extents(5)+1

		allocate(packbuf(3,ncxl,ncyl,nczl),testbuf(3,ncxl,ncyl,nczl))

		! Populate dummy packbuf
		do ixyz = 1,3
		do icell=1,extents(2)-extents(1)
		do jcell=1,extents(4)-extents(3)
		do kcell=1,extents(6)-extents(5)

			packbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                       			  1000*jcell + &
			                    			  1000000*kcell

		end do
		end do
		end do
		end do


	else if (realm .eq. cfd_realm) then	

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,extents)

		ncxl = extents(2)-extents(1)+1
		ncyl = extents(4)-extents(3)+1
		nczl = extents(6)-extents(5)+1

		allocate(packbuf(3,ncxl,ncyl,nczl),testbuf(3,ncxl,ncyl,nczl))

		! Populate dummy packbuf
		do ixyz = 1,3
		do icell=1,extents(2)-extents(1)
		do jcell=1,extents(4)-extents(3)
		do kcell=1,extents(6)-extents(5)

			packbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                       			  1000*jcell + &
			                    			  1000000*kcell

		end do
		end do
		end do
		end do

	end if

	!print*,'Test pack',rank_world,myid_cart,realm,olap_mask(rank_world), extents, coord

	call CPL_pack(packbuf,outbuf,realm)
	call CPL_unpack(outbuf,testbuf,realm)

	print'(a,3f10.5)', 'Error in pack/unpack = ', maxval(testbuf-packbuf), & 
												  minval(testbuf-packbuf), & 
												     sum(testbuf-packbuf)

end subroutine test_packing

! ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬ஜ۩۞۩ஜ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
! Test the send and recv routines from coupler

subroutine test_send_recv_MD2CFD
	use CPL 
	implicit none

	logical	:: send_flag,recv_flag
	integer :: ncxl,ncyl,nczl,ixyz,icell,jcell,kcell
	integer	:: jcmin_send,jcmax_send,jcmin_recv,jcmax_recv,npercell,coord(3),extents(6)
	double precision,dimension(:,:,:,:),allocatable	:: sendbuf,recvbuf

	npercell = 3
	jcmax_send=1; jcmin_send=1; 
	jcmax_recv = jcmax_send
	jcmin_recv = jcmin_send
	if (olap_mask(rank_world) .eqv. .false.) return

	! Test Sending from MD to CFD							   
	if (realm .eq. md_realm) then	

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_olap_extents(coord,realm,extents)

		allocate(sendbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		!print'(2a,11i7)', 'sent size',realm_name(realm),extents,size(sendbuf),shape(sendbuf)

		! Populate dummy gatherbuf
		sendbuf = -333.d0 ! 0.d0
		do ixyz = 1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)
			sendbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                                            1000*jcell + &
			                                         1000000*kcell
		end do
		end do
		end do
		end do

		call CPL_send(sendbuf,jcmax_send=jcmax_send,jcmin_send=jcmin_send,send_flag=send_flag)	

		if (send_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_send,jcmax_send
			do icell=extents(1),extents(2)
			do ixyz =1,npercell
				write(4000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				      'send MD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
				       sendbuf(ixyz,icell,jcell,kcell)
			end do
			end do
			end do
			end do
		endif

	else if (realm .eq. cfd_realm) then	 

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_proc_extents(coord,realm,extents)
		!print'(2a,8i7)', 'proc extents', realm_name(realm),rank_world,rank_cart,extents
		call CPL_olap_extents(coord,realm,extents)
		!print'(2a,8i7)', 'olap extents', realm_name(realm),rank_world,rank_cart,extents

		allocate(recvbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		!print'(2a,11i7)', 'recv size', realm_name(realm),extents,size(recvbuf),shape(recvbuf)
		recvbuf = -444.d0
		call CPL_recv(recvbuf,jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)

		if (recv_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_recv,jcmax_recv  !extents(3),extents(4)
			do icell=extents(1),extents(2)
			do ixyz =1,npercell
				!if ( recvbuf(ixyz,icell,jcell,kcell) .ne. -444.d0) then
				!print'(a,i4,a,i4,a,i4,a,i4,a,f20.1)',   &
				!      'recv CFD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
				!       recvbuf(ixyz,icell,jcell,kcell)
					write(5000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					      'recv CFD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
					       recvbuf(ixyz,icell,jcell,kcell)
				!endif

			end do
			end do
			end do
			end do
		endif
	end if								   

	! if (realm .eq.  md_realm) write(4000+myid_world,*),myid_world, 'BUF=', sendbuf
	! if (realm .eq. cfd_realm) write(5000+myid_world,*),myid_world, 'BUF=', recvbuf
	
end subroutine test_send_recv_MD2CFD


! ۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩
! Test Sending from MD to CFD

subroutine test_send_recv_CFD2MD
	use CPL
	implicit none

	logical	:: send_flag,recv_flag
	integer	:: jcmin_send,jcmax_send,jcmin_recv,jcmax_recv
	integer :: ncxl,ncyl,nczl,ixyz,icell,jcell,kcell,npercell,coord(3),extents(6)
	double precision,dimension(:,:,:,:),allocatable	:: sendbuf,recvbuf

	npercell = 3
	jcmax_send=1; jcmin_send=1; 
	jcmax_recv = jcmax_send
	jcmin_recv = jcmin_send
	if (olap_mask(rank_world) .eqv. .false.) return

	! Test Sending from CFD to MD							   
	if (realm .eq. md_realm) then		   

		coord = (/iblock_realm,jblock_realm,kblock_realm /)
		call CPL_olap_extents(coord,realm,extents)

		allocate(recvbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))
		recvbuf = -444

		!print*, 'recv size', realm_name(realm),extents, size(recvbuf),shape(recvbuf)
		call CPL_recv(recvbuf,jcmax_recv=1,jcmin_recv=1,recv_flag=recv_flag)   

		if (recv_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_send,jcmax_send
			do icell=extents(1),extents(2)
			do ixyz = 1,npercell
				write(2000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
			      	'recv MD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
			      	 recvbuf(ixyz,icell,jcell,kcell)
			end do
			end do
			end do
			end do
		endif

	else if (realm .eq. cfd_realm) then	   

		coord = (/iblock_realm,jblock_realm,kblock_realm /)
		call CPL_olap_extents(coord,realm,extents)
		allocate(sendbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		do ixyz =1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)
			sendbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*(icell) + &
			                       			  1000*(jcell) + &
			                    			  1000000*(kcell)

		end do
		end do
		end do
		end do

		!print*, 'sent size',realm_name(realm),3*ncxl*ncyl*nczl,size(sendbuf)
		call CPL_send(sendbuf,jcmax_send=1,jcmin_send=1,send_flag=send_flag)

		do kcell=extents(5),extents(6)
		do jcell=jcmin_send,jcmax_send
		do icell=extents(1),extents(2)
		do ixyz = 1,npercell
			write(9000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
			      'send CFD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
			       sendbuf(ixyz,icell,jcell,kcell)
		end do
		end do
		end do
		end do
	end if								   
	
end subroutine test_send_recv_CFD2MD

! ++ UNINTERESTING ++ ========================================================
subroutine initialise
	use mpi
	implicit none

	integer :: ierr
	
	call MPI_init(ierr)

end subroutine initialise

subroutine finalise
	use mpi
	implicit none
	
	integer :: ierr
	call MPI_finalize(ierr)

end subroutine finalise

subroutine barrier
	use mpi
	implicit none
	
	integer :: ierr	
	call MPI_barrier(MPI_COMM_WORLD,ierr)

end subroutine barrier

subroutine lasterrorcheck
	use mpi 
	implicit none
	
	integer :: ierr
	integer :: resultlen
	character*12 err_buffer

	call MPI_Error_string(ierr,err_buffer,resultlen,ierr)
	print*, err_buffer

end subroutine lasterrorcheck
