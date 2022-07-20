subroutine save_pdfi_to_sdf(output,input,sharpnum,as,bs,cs,ds,msp1,nsp1,index,index_total,tjul0,tjul1)

    !Outputs sdf files containing the CGEM maps collected through python
    !To be called iteratively - sequentially appends each timestep 

    integer, intent(in) :: msp1,nsp1,index,index_total,sharpnum
    real*8, intent(in) :: as,bs,cs,ds,tjul0(index_total),tjul1(index_total)
    double precision, intent(in) :: input(nsp1,msp1,18)
    double precision, intent(out) :: output(400,400,3)


    integer*8 :: dims(20)
    integer :: nlonpad,nlatpad,mpadb,mpadt,npadl,npadr,ms,mi,ns,ni,m,n


    real*8, allocatable :: blon0(:,:),blon1(:,:)
    real*8, allocatable :: blat0(:,:),blat1(:,:)
    real*8, allocatable :: brll0(:,:),brll1(:,:)
    real*8, allocatable :: elonpdfi(:,:),elatpdfi(:,:)
    real*8, allocatable :: delondr(:,:),delatdr(:,:)
    real*8, allocatable :: erllpdfi(:,:),erllind(:,:)
    real*8, allocatable :: srll(:,:),hmll(:,:)
    real*8, allocatable :: mcoell(:,:),mcoll(:,:),mcell(:,:)
    real*8, allocatable :: mtell(:,:),mpell(:,:)
    real*8, allocatable :: brtp(:,:),d3(:,:,:),d3i(:,:,:),eth(:,:),eph(:,:)

    real*8 :: a,b,c,d,tjulhalf
    real*8 :: srtot, hmtot
    real*8 :: times(index_total)
    real*8 :: bthr, rsun
    character*80 :: target_variable
    character*80 :: outfile_times
    character*80 :: outfile_all
    character*2 :: index_of_file
    character*4 :: num
    
    logical :: rotate=.true. !.true. sets orientation to the same as for chuns prototype


    pi=4.0d0*atan(1.0d0)
    bthr=200.0d0
    rsun=6.95d5

    write(num,'(I4.4)') sharpnum

    ! Location and name of output file contianing the times of each snapshot
    outfile_times='test_output_times_'//num//'.sdf'

    ! Location and name of output file that will contain the individual snapshots
    outfile_all='test_output_all_'//num//'.sdf'
    
    ! Convert integer to string
    write(index_of_file,'(I2.2)') index+1
    

    ! Remove already existing file name
    if (index==0) then 
        call sdf_rm_f77(outfile_times)
        call sdf_rm_f77(outfile_all)
    end if

          
    !Initial guess for padding around E field maps   
    nlatpad=50
    nlonpad=50

    !Desired resolution of maps to feed into amrvac
    mi=384
    ni=384

    allocate(blon0(ni+1,mi),blon1(ni+1,mi))
    allocate(blat0(ni,mi+1),blat1(ni,mi+1))
    allocate(brll0(ni,mi),brll1(ni,mi))
    allocate(elonpdfi(ni,mi+1))
    allocate(elatpdfi(ni+1,mi))
    allocate(delondr(ni,mi+1))
    allocate(delatdr(ni+1,mi))
    allocate(erllpdfi(ni+1,mi+1))
    allocate(erllind(ni+1,mi+1))
    allocate(srll(ni,mi),hmll(ni,mi))
    allocate(mcoell(ni+1,mi+1),mcoll(ni-1,mi-1),mcell(ni,mi))
    allocate(mtell(ni,mi+1),mpell(ni+1,mi))
    allocate(brtp(mi,ni))
    allocate(d3i(ni+1,mi+1,18))
    allocate(eth(mi,ni+1),eph(mi+1,ni))

    !Maps fed in are assumed to be on a staggered grid
    ns=nsp1-1
    ms=msp1-1

    !deduce padding based on initial guess
    call pad_int_gen_ss(ms,ns,nlatpad,nlonpad,mpadb,mpadt,npadl,npadr)
    !Apply padding to m,n,a,b,c,d variables
    call pad_abcd_as_ss(ms,ns,mpadb,mpadt,npadl,npadr,as,bs,cs,ds,m,n,a,b,c,d)

    !Print statements for visual check of padding procedure
    write(6,*) 'before pad_abcd_as_ss call: m,n = ',ms,ns
    write(6,*) 'after pad_abcd_as_ss call: m,n = ',m,n
    write(6,*) 'before pad_abcd_as_ss call: a,b,c,d = ',as,bs,cs,ds
    write(6,*) 'after pad_abcd_as_ss call: a,b,c,d = ',a,b,c,d
    if (index==0) then
        write(6,*) 'xprobmin1 = ', a/pi*0.5d0
        write(6,*) 'xprobmax1 = ', b/pi*0.5d0
        write(6,*) 'xprobmin2 == ', 0.0d0
        write(6,*) 'xprobmax2 = ', abs(c/pi*0.5d0)+abs(d/pi*0.5d0)
    end if
    
    allocate(d3(n+1,m+1,18))

    do id=1,18
        !Apply padding to input maps
        call add_padding_as_ss(ms,ns,mpadb,mpadt,npadl,npadr,input(:,:,id),d3(:,:,id))
    end do

    !Interpolate the padded input maps to the desired resolution
    call interp_hmidata_3d_ll(m,n,d3,mi,ni,d3i)

    !Pass maps through the pdfi wrapper to calculate E field inversion
    call pdfi_wrapper4jsoc_ss(mi,ni,rsun,a,b,c,d,bthr,d3i(:,:,1),&
        d3i(:,:,2),d3i(:,:,3),d3i(:,:,4),d3i(:,:,5),d3i(:,:,6),&
        d3i(:,:,7),d3i(:,:,8),d3i(:,:,9),d3i(:,:,10),d3i(:,:,11),&
        d3i(:,:,12),d3i(:,:,13),d3i(:,:,14),d3i(:,:,15),d3i(:,:,16),&
        d3i(:,:,17),d3i(:,:,18),tjul0,tjul1,blon0,blat0,brll0,blon1,&
        blat1,brll1,elonpdfi,elatpdfi,erllpdfi,erllind,delondr,delatdr,&
        srll,srtot,hmll,hmtot,mcoell,mcoll,mcell,mtell,mpell,tjulhalf)

    if (rotate) then
        !store brtp in theta-phi order i.e., following amrvac convention
        do i1=1,mi
            do i2=1,ni
              brtp(i1,i2)=brll0(i2,mi+1-i1)
            end do
        end do
        dims(1)=mi
        dims(2)=ni
        target_variable='brtp'//index_of_file
        call sdf_write_f77(outfile_all,target_variable,'f',8,2,dims,brtp)
    
        call ehyeell2tp_ss(mi,ni,elonpdfi,elatpdfi,eth,eph)
    else
        brtp=brll0
        
        dims(1)=mi
        dims(2)=ni
        call sdf_write_f77(outfile_all,target_variable,'f',8,2,dims,brtp)

        eph=elonpdfi
        eth=elatpdfi
    end if

!  first, convert E (V/cm) back to cE (G cm/sec)
    eth=eth*1.d8
    eph=eph*1.d8
    dims(1)=mi
    dims(2)=ni+1
    target_variable='eth'//index_of_file
    call sdf_write_f77(outfile_all,target_variable,'f',8,2,dims,eth)
    dims(1)=mi+1
    dims(2)=ni
    target_variable='eph'//index_of_file
    call sdf_write_f77(outfile_all,target_variable,'f',8,2,dims,eph)


    if (index==0) then
        dims(1)=index_total
        times(1)=tjul0(1)*86400
        times(2:dims(1))=tjul1(1:dims(1)-1)*86400
        print*, times
        call sdf_write_f77(outfile_times,'times','f',8,1,dims,times)
    end if

    print*, "Passed Successfully - Outfile: ", outfile_all
    print*, "Iteration: ",index_of_file

    !output component maps to be checked in python
    output(1:mi,1:ni,1)=brtp
    output(1:mi,1:ni+1,2)=eth
    output(1:mi+1,1:ni,3)=eph    


end subroutine
