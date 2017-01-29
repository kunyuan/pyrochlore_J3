SUBROUTINE def_lattice
    implicit none
    integer :: x, y, z, k, xm, ym, zm, xp, yp, zp, tmp,sit
    integer :: ix, iy, iz, ixp, iyp, izp, ixm, iym, izm, vec(3), site(4), itmp

    !!!!************DEFINITION FOR LATTICES****************************
    if(LatticeName=='Cubic') then
      	Sub(:) = 1
        if(Dim==2) then
        !RealVectors
        LatticeVector(1,:) = (/1.0, 0.0/)
        LatticeVector(2,:) = (/0.0, 1.0/)
        SubVector(1,:) = (/0.0, 0.0/)

        sit=0
        do iy=1,L(2)
          do ix=1,L(1)
            RealVector(sit,:)=(ix-1)*LatticeVector(1,:)+(iy-1)*LatticeVector(2,:)
            sit = sit+1
          enddo
        enddo

	  k = 0
	  do y = 1,L(2)
	    do x = 1,L(1)
	      k = k + 1

	      xm =  -1 ;        if(x== 1) xm = L(1)-1
	      xp =  +1 ;        if(x==L(1)) xp = 1-L(1)

	      ym =  -L(1);      if(y== 1) ym = Vol-L(1)
	      yp =  +L(1);      if(y==L(2)) yp = L(1)-Vol

	      Ngs(1,k) = k+xp;    Ngs(4,k) = k+xm
	      Ngs(2,k) = k+yp;    Ngs(3,k) = k+ym

		    !            2
		    !            |
		    !       4---   ---1
		    !            |
		    !            3
	    enddo
	  enddo 

	  Back(1, 1) = 4; Back(4, 1) = 1;
	  Back(2, 1) = 3; Back(3, 1) = 2;

	  !-- auxillary variables to measure wrapping probability-----------
	  dr(:,:,:)=0
	  dr(1,1, 1) = 1;  dr(4,1, 1) =-1
	  dr(2,1, 2) = 1;  dr(3,1, 2) =-1
	else if(Dim==3) then
	  !RealVectors
	  LatticeVector(1,:) = (/1.0, 0.0, 0.0/)
	  LatticeVector(2,:) = (/0.0, 1.0, 0.0/)
	  LatticeVector(3,:) = (/0.0, 0.0, 1.0/)
	  SubVector(1,:) = (/0.0, 0.0, 0.0/)
	  sit=0
	  do iz=1, L(3)
	    do iy=1,L(2)
	      do ix=1,L(1)
	        RealVector(sit,:)=(ix-1)*LatticeVector(1,:)+(iy-1)*LatticeVector(2,:)+(iz-1)*LatticeVector(3,:)
	        sit = sit+1
	      enddo
	    enddo
	  enddo

	  k = 0
	  do z = 1,L(3)
	    do y = 1,L(2)
	      do x = 1,L(1)
		k = k + 1

		xm =  -1 ;      if(x== 1) xm = L(1)-1
		xp =  +1 ;      if(x==L(1)) xp = 1-L(1)

		ym =  -L(1);      if(y== 1) ym = L(1)*L(2)-L(1)
		yp =  +L(1);      if(y==L(2)) yp = L(1)-L(1)*L(2)

		zm =  -L(1)*L(2);      if(z== 1) zm = Vol-L(1)*L(2)
		zp =  +L(1)*L(2);      if(z==L(3)) zp = L(1)*L(2)-Vol

		Ngs(1,k) = k+xp;    Ngs(4,k) = k+xm
		Ngs(2,k) = k+yp;    Ngs(5,k) = k+ym
		Ngs(3,k) = k+zp;    Ngs(6,k) = k+zm
	    enddo
       	  enddo 
	enddo

        Back(1, 1) = 4; Back(4, 1) = 1;
        Back(2, 1) = 5; Back(5, 1) = 2;
        Back(3, 1) = 6; Back(6, 1) = 3;

	!-- auxillary variables to measure wrapping probability-----------
	dr(:,:,:) = 0
	dr(1,1, 1) = 1;   dr(4,1, 1) = -1
	dr(2,1, 2) = 1;   dr(5,1, 2) = -1
	dr(3,1, 3) = 1;   dr(6,1, 3) = -1
      endif

    else if(LatticeName=='Pyrochlore') then
      !RealVectors
      LatticeVector(1,:) = (/0.0, 1.0, 1.0/)
      LatticeVector(2,:) = (/1.0, 0.0, 1.0/)
      LatticeVector(3,:) = (/1.0, 1.0, 0.0/)
      SubVector(1,:) = (/0.0, 0.0, 0.0/)
      SubVector(2,:) = (/0.0, 0.5, 0.5/)
      SubVector(3,:) = (/0.5, 0.0, 0.5/)
      SubVector(4,:) = (/0.5, 0.5, 0.0/)

      do i = 1, 4
	site(i) = i
      enddo
      do iz=1,L(3)
        do iy=1,L(2)
       	  do ix=1,L(1)
	    RealVector(site(1),:)=(ix-1)*LatticeVector(1,:)+(iy-1)*LatticeVector(2,:)+(iz-1)*LatticeVector(3,:)
	    RealVector(site(2),:)=RealVector(site(1),:)+SubVector(2,:)
	    RealVector(site(3),:)=RealVector(site(1),:)+SubVector(3,:)
	    RealVector(site(4),:)=RealVector(site(1),:)+SubVector(4,:)
	    do i = 1, 4
		site(i) = site(i)+4
	    enddo
	  enddo
        enddo
      enddo

      !directions
      Back(1, 1) = 4; Back(4, 2) = 1;
      Back(2, 1) = 4; Back(4, 3) = 2;
      Back(3, 1) = 4; Back(4, 4) = 3;
      Back(4, 1) = 1; Back(1, 2) = 4;
      Back(5, 1) = 1; Back(1, 3) = 5;
      Back(6, 1) = 1; Back(1, 4) = 6;

      Back(2, 2) = 5; Back(5, 3) = 2;
      Back(3, 2) = 5; Back(5, 4) = 3;

      Back(5, 2) = 2; Back(2, 3) = 5;
      Back(6, 2) = 2; Back(2, 4) = 6;

      Back(3, 3) = 6; Back(6, 4) = 3;
      Back(6, 3) = 3; Back(3, 4) = 6;

      !lattice connections
      do i = 1, 4
        site(i) = i
	Sub(site(i)) = i
      enddo

      do iz=1,L(3)
	do iy=1,L(2)
	  do ix=1,L(1)
            ixp = ix+1
	    ixm = ix-1
	    iyp = iy+1
	    iym = iy-1
	    izp = iz+1
	    izm = iz-1

            if(ixp>L(1)) ixp = ixp-L(1)
	    if(ixm==0) ixm = L(1)
	    if(iyp>L(2)) iyp = iyp-L(2)
	    if(iym==0) iym = L(2)
	    if(izp>L(3)) izp = izp-L(3)
	    if(izm==0) izm = L(3)

	    !!!!!NEAREST NEIGHBOR!!!!!!!!!!!!!!!!!!!!!!
	    do i = 1, 3
	      Ngs(i, site(1)) = site(i+1)
	      Ngs(back(i, 1),site(i+1)) = site(1)
	    enddo

    	    call GetSite(itmp, (/ixp, iy, iz/), 1)
	    Ngs(1, site(2))  = itmp
	    Ngs(back(1, 2), itmp) = site(2)

     	    do i = 2, 3
	      Ngs(i, site(2)) = site(i+1)
	      Ngs(back(i, 2), site(i+1)) = site(2)
	    enddo

    	    call GetSite(itmp, (/ix, iyp, iz/), 1)
	    Ngs(1, site(3)) =  itmp
	    Ngs(back(1, 3), itmp) = site(3)

    	    call GetSite(itmp, (/ixm, iyp, iz/), 2)
	    Ngs(2, site(3)) = itmp
	    Ngs(back(2, 3), itmp) = site(3)

	    Ngs(3, site(3)) = site(4)
	    Ngs(back(3, 3), site(4)) = site(3)

    	    call GetSite(itmp, (/ix, iy, izp/), 1)
	    Ngs(1, site(4)) = itmp
	    Ngs(back(1, 4), itmp) = site(4)

    	    call GetSite(itmp, (/ixm, iy, izp/), 2)
	    Ngs(2, site(4)) = itmp
	    Ngs(back(2, 4), itmp) = site(4)

	    call GetSite(itmp, (/ix, iym, izp/), 3)
	    Ngs(3, site(4)) = itmp
	    Ngs(back(3, 4), itmp) = site(4)

	    do i = 1, NSub
		site(i) = site(i)+4
		Sub(site(i)) = i
	    enddo
	  enddo
	enddo
      enddo

      !!!!!Third NEAREST NEIGHBORS, connected!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, 4
	  Back(7, i) = 10;   Back(8, i) = 11;  Back(9, i) = 12
	  Back(10,i) = 7;   Back(11, i) = 8;   Back(12,i) = 9
      enddo

      do i = 1, 4
        site(i) = i
      enddo

      do iz=1,L(3)
	do iy=1,L(2)
	  do ix=1,L(1)
            ixp = ix+1
	    ixm = ix-1
	    iyp = iy+1
	    iym = iy-1
	    izp = iz+1
	    izm = iz-1

            if(ixp>L(1)) ixp = ixp-L(1)
	    if(ixm==0) ixm = L(1)
	    if(iyp>L(2)) iyp = iyp-L(2)
	    if(iym==0) iym = L(2)
	    if(izp>L(3)) izp = izp-L(3)
	    if(izm==0) izm = L(3)

    	    call GetSite(itmp, (/ixp, iy, iz/), 1);   Ngs(7, site(1))= itmp
    	    call GetSite(itmp, (/ix, iyp, iz/), 1);   Ngs(8, site(1))= itmp
    	    call GetSite(itmp, (/ix, iy, izp/), 1);   Ngs(9, site(1))= itmp
    	    call GetSite(itmp, (/ixm, iy, iz/), 1);   Ngs(10, site(1))= itmp
    	    call GetSite(itmp, (/ix, iym, iz/), 1);   Ngs(11, site(1))= itmp
    	    call GetSite(itmp, (/ix, iy, izm/), 1);   Ngs(12, site(1))= itmp

    	    call GetSite(itmp, (/ixm, iy, iz/), 2);   Ngs(7, site(2))= itmp
    	    call GetSite(itmp, (/ixm, iyp, iz/), 2);  Ngs(8, site(2))= itmp
    	    call GetSite(itmp, (/ixm, iy, izp/),2);   Ngs(9, site(2))= itmp
    	    call GetSite(itmp, (/ixp, iy, iz/),2);    Ngs(10, site(2))= itmp
    	    call GetSite(itmp, (/ixp, iym, iz/),2);   Ngs(11, site(2))= itmp
    	    call GetSite(itmp, (/ixp, iy, izm/),2);   Ngs(12, site(2))= itmp

    	    call GetSite(itmp, (/ix, iym, iz/),  3);   Ngs(7, site(3))= itmp
    	    call GetSite(itmp, (/ixp, iym, iz/), 3);   Ngs(8, site(3))= itmp
    	    call GetSite(itmp, (/ix, iym, izp/), 3);   Ngs(9, site(3))= itmp
    	    call GetSite(itmp, (/ix, iyp, iz/),  3);   Ngs(10, site(3))= itmp
    	    call GetSite(itmp, (/ixm, iyp, iz/), 3);   Ngs(11, site(3))= itmp
    	    call GetSite(itmp, (/ix, iyp, izm/), 3);   Ngs(12, site(3))= itmp

    	    call GetSite(itmp, (/ix, iy,  izm/), 4);   Ngs(7, site(4))= itmp
    	    call GetSite(itmp, (/ixp, iy, izm/), 4);   Ngs(8, site(4))= itmp
    	    call GetSite(itmp, (/ix, iyp, izm/), 4);   Ngs(9, site(4))= itmp
    	    call GetSite(itmp, (/ix, iy,  izp/), 4);   Ngs(10, site(4))= itmp
    	    call GetSite(itmp, (/ixm, iy, izp/), 4);   Ngs(11, site(4))= itmp
    	    call GetSite(itmp, (/ix, iym, izp/), 4);   Ngs(12, site(4))= itmp

	    do i = 1, NSub
		site(i) = site(i)+4
	    enddo
	  enddo
        enddo
      enddo

      !!!!!Second NEAREST NEIGHBORS, connected!!!!!!!!!!!!!!!!!!!!!!
      !do i = 1, 4
        !site(i) = i
      !enddo

      !do iz=1,L(3)
	!do iy=1,L(2)
	  !do ix=1,L(1)
            !ixp = ix+1
	    !ixm = ix-1
	    !iyp = iy+1
	    !iym = iy-1
	    !izp = iz+1
	    !izm = iz-1

            !if(ixp>L(1)) ixp = ixp-L(1)
	    !if(ixm==0) ixm = L(1)
	    !if(iyp>L(2)) iyp = iyp-L(2)
	    !if(iym==0) iym = L(2)
	    !if(izp>L(3)) izp = izp-L(3)
	    !if(izm==0) izm = L(3)

                !call GetSite(itmp, (/ix, iym, iz/), 2);   Ngs(7, site(1))= itmp
                !call GetSite(itmp, (/ix, iy, izm/), 2);   Ngs(8, site(1))= itmp
                !call GetSite(itmp, (/ixm, iyp, iz/), 2);  Ngs(9, site(1))= itmp
                !call GetSite(itmp, (/ixm, iy, izp/), 2);  Ngs(10, site(1))= itmp

                !call GetSite(itmp, (/ixp, iym, iz/), 1);   Ngs(7, site(2))= itmp
                !call GetSite(itmp, (/ixp, iy, izm/), 1);   Ngs(8, site(2))= itmp
                !call GetSite(itmp, (/ix, iy, izp/), 1);  Ngs(9, site(2))= itmp
                !call GetSite(itmp, (/ix, iyp, iz/), 1);  Ngs(10, site(2))= itmp
	    !Back(7, 1) = 10;   Back(8, 1) = 9;  Back(9, 1) = 7; Back(10, 1) = 8
	    !Back(10, 2) = 7;   Back(9, 2) = 8;  Back(7, 2) = 9; Back(8, 2) = 10


                !call GetSite(itmp, (/ixm, iy, iz/), 3);   Ngs(11, site(1))= itmp
                !call GetSite(itmp, (/ix, iy, izm/), 3);   Ngs(12, site(1))= itmp
                !call GetSite(itmp, (/ixp, iym, iz/), 3);  Ngs(13, site(1))= itmp
                !call GetSite(itmp, (/ix, iym, izp/), 3);  Ngs(14, site(1))= itmp

                !call GetSite(itmp, (/ixm, iyp, iz/), 1);   Ngs(7, site(3))= itmp
                !call GetSite(itmp, (/ix, iyp, izm/), 1);   Ngs(8, site(3))= itmp
                !call GetSite(itmp, (/ixp, iy, iz/), 1);  Ngs(9, site(3))= itmp
                !call GetSite(itmp, (/ix, iy, izp/), 1);  Ngs(10, site(3))= itmp
	    !Back(11, 1) = 9;   Back(12, 1) = 10;  Back(13, 1) = 7; Back(14, 1) = 8
	    !Back(9, 3) = 11;   Back(10, 3) = 12;  Back(7, 3) = 13; Back(8, 3) = 14

                !call GetSite(itmp, (/ixm, iy, iz/), 4);   Ngs(15, site(1))= itmp
                !call GetSite(itmp, (/ix, iym, iz/), 4);   Ngs(16, site(1))= itmp
                !call GetSite(itmp, (/ixp, iy, izm/), 4);  Ngs(17, site(1))= itmp
                !call GetSite(itmp, (/ix, iyp, izm/), 4);  Ngs(18, site(1))= itmp

                !call GetSite(itmp, (/ixm, iy, izp/), 1);   Ngs(7, site(4))= itmp
                !call GetSite(itmp, (/ix, iym, izp/), 1);   Ngs(8, site(4))= itmp
                !call GetSite(itmp, (/ixp, iy, iz/), 1);    Ngs(9, site(4))= itmp
                !call GetSite(itmp, (/ix, iyp, iz/), 1);    Ngs(10, site(4))= itmp
	    !Back(15, 1) = 9;   Back(16, 1) = 10;  Back(17, 1) = 7; Back(18, 1) = 8
	    !Back(9, 4) = 15;   Back(10, 4) = 16;  Back(7, 4) = 17; Back(8, 4) = 18

                !call GetSite(itmp, (/ixp, iy, iz/), 3);   Ngs(11, site(2))= itmp
                !call GetSite(itmp, (/ixp, iy, izm/), 3);   Ngs(12, site(2))= itmp
                !call GetSite(itmp, (/ix, iym, iz/), 3);  Ngs(13, site(2))= itmp
                !call GetSite(itmp, (/ix, iym, izp/), 3);  Ngs(14, site(2))= itmp

                !call GetSite(itmp, (/ix, iyp, iz/), 2);   Ngs(11, site(3))= itmp
                !call GetSite(itmp, (/ix, iyp, izm/), 2);   Ngs(12, site(3))= itmp
                !call GetSite(itmp, (/ixm, iy, iz/), 2);  Ngs(13, site(3))= itmp
                !call GetSite(itmp, (/ixm, iy, izp/), 2);  Ngs(14, site(3))= itmp
	    !Back(11, 2) = 13;   Back(12, 2) = 14;  Back(13, 2) = 11; Back(14, 2) = 12
	    !Back(13, 3) = 11;   Back(14, 3) = 12;  Back(11, 3) = 13; Back(12, 3) = 14

                !call GetSite(itmp, (/ixp, iym, iz/), 4);   Ngs(15, site(2))= itmp
                !call GetSite(itmp, (/ixp, iy, iz/), 4);   Ngs(16, site(2))= itmp
                !call GetSite(itmp, (/ix, iy, izm/), 4);  Ngs(17, site(2))= itmp
                !call GetSite(itmp, (/ix, iyp, izm/), 4);  Ngs(18, site(2))= itmp

                !call GetSite(itmp, (/ix, iy, izp/), 2);    Ngs(11, site(4))= itmp
                !call GetSite(itmp, (/ix, iym, izp/), 2);   Ngs(12, site(4))= itmp
                !call GetSite(itmp, (/ixm, iy, iz/), 2);    Ngs(13, site(4))= itmp
                !call GetSite(itmp, (/ixm, iyp, iz/), 2);   Ngs(14, site(4))= itmp
	    !Back(15, 2) = 14;   Back(16, 2) = 13;  Back(17, 2) = 11; Back(18, 2) = 12
	    !Back(14, 4) = 15;   Back(13, 4) = 16;  Back(11, 4) = 17; Back(12, 4) = 18


                !call GetSite(itmp, (/ix, iyp, iz/), 4);   Ngs(15, site(3))= itmp
                !call GetSite(itmp, (/ixm, iyp, iz/), 4);   Ngs(16, site(3))= itmp
                !call GetSite(itmp, (/ix, iy, izm/), 4);  Ngs(17, site(3))= itmp
                !call GetSite(itmp, (/ixp, iy, izm/), 4);  Ngs(18, site(3))= itmp

                !call GetSite(itmp, (/ix, iy, izp/), 3);    Ngs(15, site(4))= itmp
                !call GetSite(itmp, (/ixm, iy, izp/), 3);   Ngs(16, site(4))= itmp
                !call GetSite(itmp, (/ix, iym, iz/), 3);    Ngs(17, site(4))= itmp
                !call GetSite(itmp, (/ixp, iym, iz/), 3);   Ngs(18, site(4))= itmp
	    !Back(15, 3) = 17;   Back(16, 3) = 18;  Back(17, 3) = 15; Back(18, 3) = 16
	    !Back(17, 4) = 15;   Back(18, 4) = 16;  Back(15, 4) = 17; Back(16, 4) = 18


	    !do i = 1, NSub
		!site(i) = site(i)+4
	    !enddo
	  !enddo
        !enddo
      !enddo

      !-- auxillary variables to measure wrapping probability-----------
      !-- dr(nnb, sublattice of the target site, x/y/z)-----------------
      dr(:,:,:) = 0
      dr(1,1, 1) = 1;   dr(1,1, 2) =  0;   dr(1,1, 3) =  0
      dr(4,1, 1) =-1;   dr(4,1, 2) =  0;   dr(4,1, 3) =  0 
      dr(2,1, 1) = 0;   dr(2,1, 2) =  1;   dr(2,1, 3) =  0
      dr(5,1, 1) = 0;   dr(5,1, 2) = -1;   dr(5,1, 3) =  0 
      dr(3,1, 1) = 0;   dr(3,1, 2) =  0;   dr(3,1, 3) =  1    
      dr(6,1, 1) = 0;   dr(6,1, 2) =  0;   dr(6,1, 3) = -1     

      dr(1,2, 1) = 1;   dr(1,2, 2) =  0;   dr(1,2, 3) =  0
      dr(4,2, 1) =-1;   dr(4,2, 2) =  0;   dr(4,2, 3) =  0 
      dr(2,2, 1) =-1;   dr(2,2, 2) =  1;   dr(2,2, 3) =  0
      dr(5,2, 1) = 1;   dr(5,2, 2) = -1;   dr(5,2, 3) =  0 
      dr(3,2, 1) =-1;   dr(3,2, 2) =  0;   dr(3,2, 3) =  1
      dr(6,2, 1) = 1;   dr(6,2, 2) =  0;   dr(6,2, 3) = -1 

      dr(1,3, 1) = 0;   dr(1,3, 2) =  1;   dr(1,3, 3) =  0
      dr(4,3, 1) = 0;   dr(4,3, 2) = -1;   dr(4,3, 3) =  0 
      dr(2,3, 1) =-1;   dr(2,3, 2) =  1;   dr(2,3, 3) =  0
      dr(5,3, 1) = 1;   dr(5,3, 2) = -1;   dr(5,3, 3) =  0 
      dr(3,3, 1) = 0;   dr(3,3, 2) = -1;   dr(3,3, 3) =  1
      dr(6,3, 1) = 0;   dr(6,3, 2) =  1;   dr(6,3, 3) = -1 

      dr(1,4, 1) = 0;   dr(1,4, 2) =  0;   dr(1,4, 3) =  1
      dr(4,4, 1) = 0;   dr(4,4, 2) =  0;   dr(4,4, 3) = -1 
      dr(2,4, 1) =-1;   dr(2,4, 2) =  0;   dr(2,4, 3) =  1
      dr(5,4, 1) = 1;   dr(5,4, 2) =  0;   dr(5,4, 3) = -1 
      dr(3,4, 1) = 0;   dr(3,4, 2) = -1;   dr(3,4, 3) =  1
      dr(6,4, 1) = 0;   dr(6,4, 2) =  1;   dr(6,4, 3) = -1 

      call PrintLattice
    else if(LatticeName=='Triangular') then
	Sub(:) = 1
        !RealVectors
        LatticeVector(1,:) = (/1.0, 0.0/)
        LatticeVector(2,:) = (/0.5, sqrt(3.0)/2.0/)
        SubVector(1,:) = (/0.0, 0.0/)

        sit=0
        do iy=1,L(2)
	    do ix=1,L(1)
	      RealVector(sit,:)=(ix-1)*LatticeVector(1,:)+(iy-1)*LatticeVector(2,:)
	      sit = sit+1
	    enddo
	enddo

        k = 0
        do y = 1,L(2)
	    do x = 1,L(1)
	      k = k + 1

	      xm =  -1 ;        if(x== 1) xm = L(1)-1
	      xp =  +1 ;        if(x==L(1)) xp = 1-L(1)

	      ym =  -L(1);      if(y== 1) ym = Vol-L(1)
	      yp =  +L(1);      if(y==L(2)) yp = L(1)-Vol

	      Ngs(1,k) = k+xp;       Ngs(6,k) = k+xm
	      Ngs(2,k) = k+yp;       Ngs(5,k) = k+ym
	      Ngs(3,k) = k+yp+xm;    Ngs(4,k) = k+ym+xp

	enddo
      enddo 

      Back(1, 1) = 6;   Back(6, 1) = 1 
      Back(2, 1) = 5;   Back(5, 1) = 2
      Back(3, 1) = 4;   Back(4, 1) = 3

      !-- auxillary variables to measure wrapping probability-----------
      dr(:,:,:)=0
      dr(1,1, 1) = 1;  dr(6, 1, 1) = -1
      dr(2,1, 2) = 1;  dr(5, 1, 2) = -1
      dr(3,1, 1) =-1;  dr(4, 1, 1) =  1
      dr(3,1, 2) = 1;  dr(4, 1, 2) = -1
      endif
  END SUBROUTINE def_lattice

  SUBROUTINE GetVector(Site, Vector, Sub)
    implicit none
    integer, intent(in) :: Site
    integer, intent(out) :: Vector(1:Dim), Sub
    integer :: tmp
    if(LatticeName=='Cubic') then
	if(Dim==2) then
	    Vector(2) = (Site-1)/L(1)+1
	    Vector(1) = Site-(Vector(2)-1)*L(1)
	    Sub = 1
	else if(Dim==3) then
	    Vector(3) = (Site-1)/L(1)/L(2)+1
	    Vector(2) = (Site-1-(Vector(3)-1)*L(1)*L(2))/L(1)+1
	    Vector(1) = Site-(Vector(3)-1)*L(1)*L(2)-(Vector(2)-1)*L(1)
	    Sub = 1
	endif
    else if(LatticeName=='Pyrochlore') then
	Sub = mod(Site, NSub)
	if(Sub==0) Sub = 4
	tmp = (Site-1)/4
	Vector(3) = tmp/L(1)/L(2)+1
	Vector(2) = (tmp-(Vector(3)-1)*L(1)*L(2))/L(1)+1
	Vector(1) = tmp-(Vector(3)-1)*L(1)*L(2)-(Vector(2)-1)*L(1)+1
    else if(LatticeName=='Triangular') then
	Vector(2) = (Site-1)/L(1)+1
	Vector(1) = Site-(Vector(2)-1)*L(1)
	Sub = 1
    endif
  END SUBROUTINE

  SUBROUTINE GetSite(Site, Vector, Sub)
    implicit none
    integer, intent(out) :: Site
    integer, intent(in) :: Vector(1:Dim), Sub
    if(LatticeName=='Cubic') then
	if(Dim==2) then
	    Site = Vector(1)+(Vector(2)-1)*L(1)
	else if(Dim==3) then
	    Site = (Vector(3)-1)*L(1)*L(2)+(Vector(2)-1)*L(1)+Vector(1)
	endif
    else if(LatticeName=='Pyrochlore') then
	Site = ((Vector(3)-1)*L(1)*L(2)+(Vector(2)-1)*L(1)+Vector(1)-1)*NSub+Sub
    else if(LatticeName=='Triangular') then
	Site = Vector(1)+(Vector(2)-1)*L(1)
    endif
  END SUBROUTINE

  SUBROUTINE PrintLattice()
    implicit none
    integer :: i, j, k, nb
    double precision :: deltar(3)
    if(LatticeName=="Pyrochlore") then
	open(18,file='coord.txt');  
	write(18,*) "{'Points': ["
	DO j=1, Vol 
	    write(18,*) "[(", RealVector(j,1), ",", RealVector(j,2), ",", RealVector(j,3), "), ", Sub(j), "],"
	END DO 
	write(18,*) "], "
	write(18, *) "'Interactions': ["
	DO j=1, Vol 
	    DO nb = 1,nnb
		write(18,*) "[", j, ", ",  Ngs(nb, j), "],"
		deltar(:) = RealVector(Ngs(nb, j),:)-RealVector(j, :)
		do i = 1, Dim
		    if(deltar(i)>=L(i)-0.5) then
			deltar(i) = deltar(i) -L(i)
		    endif
		    if(deltar(i)<=-(L(i)-0.5)) then
			deltar(i) = deltar(i) +L(i)
		    endif
		    deltar(i) = deltar(i)-sum(dr(nb, Sub(j), :)*LatticeVector(:, i)/2.0)
		    if(deltar(i)/=0) then
			print *, j, nb, i, deltar(i)
		    endif
		enddo
	    ENDDO
	END DO 
	write(18,*) "], "
	write(18,*) "}"
	close(18)

    else if(LatticeName=="Triangular") then
	open(18,file='coord.txt');  
	write(18,*) "{'Points': ["
	DO j=1, Vol 
	    write(18,*) "[(", RealVector(j,1), ",", RealVector(j,2), "), ", Sub(j), "],"
	END DO 
	write(18,*) "], "
	write(18, *) "'Interactions': ["
	DO j=1, Vol 
	    DO nb = 1,nnb
		write(18,*) "[", j, ", ",  Ngs(nb, j), "],"
		deltar(:) = RealVector(Ngs(nb, j),:)-RealVector(j, :)
		do i = 1, Dim
		    if(deltar(i)>=L(i)-0.5) then
			deltar(i) = deltar(i) -L(i)
		    endif
		    if(deltar(i)<=-(L(i)-0.5)) then
			deltar(i) = deltar(i) +L(i)
		    endif
		    deltar(i) = deltar(i)-sum(dr(nb, Sub(j), :)*LatticeVector(:, i)/2.0)
		    if(deltar(i)/=0) then
			print *, j, nb, i, deltar(i)
		    endif
		enddo
	    ENDDO
	END DO 
	write(18,*) "], "
	write(18,*) "}"
	close(18)
    endif
  END SUBROUTINE PrintLattice


