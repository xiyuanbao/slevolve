module constants_mod
  implicit none
  real, parameter :: pi = 3.1415926535898
  complex, parameter :: ii = (0.0,1.0)
  real, parameter :: gravConst = 6.67408E-11
end module constants_mod

module user_specs_mod
  implicit none
  integer, parameter :: str_len=200, npam=500
  integer :: unit_num=6
  character(str_len) :: inputfolder_ice='INPUT_FILES/icemodel/'
  character(str_len) :: inputfolder='INPUT_FILES/others/'
  character(str_len) :: planetfolder='INPUT_FILES/earthmodel/'
  character(str_len) :: gridfolder='INPUT_FILES/others/'
  character(str_len) :: outputfolder='OUTPUT_SLM/'
  character(str_len) :: outputfolder_ice='ICELOAD_SLM/'
  character(str_len) :: folder_coupled=''
  character(4) :: ext=''
  character(6) :: fType_in='text', fType_out='text'
  character(5) :: whichplanet='earth'
  character(str_len) :: planetmodel='prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22'
  character(str_len) :: icemodel='AISload_gd_Run85_0_0_'
  character(str_len) :: icemodel_out='iceload_out_'
  character(str_len) :: timearray='times'
  character(str_len) :: topomodel='etopo2_512_AIS_Run85_0_0'
  character(str_len) :: topo_initial='etopo2_512_AIS_Run85_0_0'
  character(str_len) :: grid_lat='GLlat_512.txt', grid_lon='GLlon_512.txt'
  integer :: norder=512, nglv=512
  real, parameter :: epsilon1=1.0E30, epsilon2=1.0E-5
  logical :: checkmarine=.false., tpw=.true., calcRG=.true.
  ! second-order load placement (build switch RAMP=1); default: step placement (Kendall et al. 2005)
  logical :: ramp=.false.
  ! axial rotation (build switch AXIAL_ROTATION=off|deg2|full): 0 = Mitrovica & Wahr (2011, eq. 15) (default),
  ! 1 = + rotational I33 (tidal k_2 acting on Lambda_20), 2 = + degree-0 trace term (elastic axial number tau_c)
  integer :: axial_rotation=0
  real :: axial_tau_c=5.360266227471394E-4
  logical :: input_times=.false., initial_topo=.true., iceVolume=.false.
  logical :: coupling=.false., patch_ice=.false.
  integer :: L_sim=2
  integer :: dt1=1,dt2=1,dt3=1,dt4=1
  integer :: Ldt1=2,Ldt2=0,Ldt3=0,Ldt4=0
end module user_specs_mod

module planets_mod
  use constants_mod
  implicit none
  real :: radius,mass,rhoi,rhow,gacc,omega,acoef,ccoef,moiA,moiC,kf
contains
  subroutine earth_init
    radius=6.371E6; mass=5.976E24; rhoi=920.0; rhow=1000.0; gacc=9.80665
    omega=7.292e-5; moiA=0.3296145*mass*radius**2; moiC=0.3307007*mass*radius**2; kf=0.9342+0.008
  end subroutine earth_init
  subroutine mars_init
    radius=3.3899E6; mass=6.4185E23; rhoi=1220.0; rhow=1000.0; gacc=3.713
    omega=7.08819118E-5; moiA=0.363914*mass*radius**2; moiC=0.365905*mass*radius**2; kf=0.899
  end subroutine mars_init
end module planets_mod

module sl_io_mod
  use user_specs_mod
  implicit none
contains
  subroutine read_sl(data_slm,filename,filepath,suffix,fext)
    real,dimension(:,:),intent(inout)::data_slm
    character(len=*),intent(in)::filename,filepath
    character(len=*),optional,intent(in)::suffix,fext
    character(len=1024)::path,xe
    xe=ext; if(present(fext)) xe=fext
    if(present(suffix)) then; path=trim(filepath)//trim(filename)//trim(suffix)//trim(xe)
    else; path=trim(filepath)//trim(filename)//trim(xe); endif
    open(unit=201,file=trim(path),form='formatted',access='sequential',status='old',action='read')
    read(201,*) data_slm
    close(201)
  end subroutine read_sl
  subroutine write_sl(data_slm,filename,filepath,suffix,fext)
    real,dimension(:,:),intent(in)::data_slm
    character(len=*),intent(in)::filename,filepath
    character(len=*),optional,intent(in)::suffix,fext
    character(len=1024)::path,xe
    xe=ext; if(present(fext)) xe=fext
    if(present(suffix)) then; path=trim(filepath)//trim(filename)//trim(suffix)//trim(xe)
    else; path=trim(filepath)//trim(filename)//trim(xe); endif
    open(unit=201,file=trim(path),form='formatted',access='sequential',status='replace',action='write')
    write(201,'(ES16.9E2)') data_slm
    close(201)
  end subroutine write_sl
  subroutine sl_get_model_res(no,ng)
    integer,intent(out)::no,ng; no=norder; ng=nglv
  end subroutine sl_get_model_res
  subroutine sl_readnl(inputfolder_ice_o,inputfolder_o,planetfolder_o,gridfolder_o,outputfolder_o,outputfolder_ice_o,folder_coupled_o,ext_o,fType_in_o,fType_out_o,planetmodel_o,icemodel_o,icemodel_out_o,timearray_o,topomodel_o,topo_initial_o,grid_lat_o,grid_lon_o,checkmarine_o,tpw_o,calcRG_o,input_times_o,initial_topo_o,iceVolume_o,coupling_o,patch_ice_o,L_sim_o,dt1_o,dt2_o,dt3_o,dt4_o,Ldt1_o,Ldt2_o,Ldt3_o,Ldt4_o,whichplanet_o)
    character(len=*),intent(out)::inputfolder_ice_o,inputfolder_o,planetfolder_o,gridfolder_o,outputfolder_o,outputfolder_ice_o,folder_coupled_o,ext_o,fType_in_o,fType_out_o,planetmodel_o,icemodel_o,icemodel_out_o,timearray_o,topomodel_o,topo_initial_o,grid_lat_o,grid_lon_o,whichplanet_o
    logical,intent(out)::checkmarine_o,tpw_o,calcRG_o,input_times_o,initial_topo_o,iceVolume_o,coupling_o,patch_ice_o
    integer,intent(out)::L_sim_o,dt1_o,dt2_o,dt3_o,dt4_o,Ldt1_o,Ldt2_o,Ldt3_o,Ldt4_o
    inputfolder_ice_o=inputfolder_ice; inputfolder_o=inputfolder; planetfolder_o=planetfolder; gridfolder_o=gridfolder
    outputfolder_o=outputfolder; outputfolder_ice_o=outputfolder_ice; folder_coupled_o=folder_coupled; ext_o=ext
    fType_in_o=fType_in; fType_out_o=fType_out; planetmodel_o=planetmodel; icemodel_o=icemodel; icemodel_out_o=icemodel_out
    timearray_o=timearray; topomodel_o=topomodel; topo_initial_o=topo_initial; grid_lat_o=grid_lat; grid_lon_o=grid_lon
    checkmarine_o=checkmarine; tpw_o=tpw; calcRG_o=calcRG; input_times_o=input_times; initial_topo_o=initial_topo
    iceVolume_o=iceVolume; coupling_o=coupling; patch_ice_o=patch_ice; L_sim_o=L_sim; dt1_o=dt1; dt2_o=dt2
    dt3_o=dt3; dt4_o=dt4; Ldt1_o=Ldt1; Ldt2_o=Ldt2; Ldt3_o=Ldt3; Ldt4_o=Ldt4; whichplanet_o=whichplanet
  end subroutine sl_readnl
end module sl_io_mod
