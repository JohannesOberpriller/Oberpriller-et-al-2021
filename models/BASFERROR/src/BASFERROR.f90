Subroutine BASFERROR(RESTART, FORTYPE, STATESPACE, RANDERR, BIAS, &
	PARAMS,MATRIX_WEATHER,CALENDAR_FERT,CALENDAR_NDEP, &
	CALENDAR_PRUNT,CALENDAR_THINT, NDAYS, NOUT, STATEVARS, &
	STATEERS, PROCERR, y)
!=====================================================
! This is the BASFOR model.
! Authors: Marcel van Oijen, mvano@ceh.ac.uk
!          David Cameron   , dcam@ceh.ac.uk
! MODEL VERSION: CONIFEROUS & DECIDUOUS
! Date: 2015-04-21
!=====================================================

use parameters
use environment
use tree
use belowgroundres
use soil
use management
use errors
implicit none

! As long as the total number of parameters stays below 100, the next line need not be changed
integer                   :: RESTART    ! Defining restart or model or not (0 : new model, 1: restart)
integer                   :: FORTYPE
integer                   :: STATESPACE ! Defining if state-space approach is used (0: no, 1: yes)
integer                   :: RANDERR ! Definin if there is a statistical error (0 : no error, 1 : error)
integer                   :: BIAS   ! Defining if there is a bias in the process or not (0: no bias, 1: bias)
integer, parameter        :: NPAR     = 100
real                      :: STATEERS(3)
real                      :: PROCERR(5)
real                      :: PARAMS(NPAR)
real                      :: STATEVARS(14)
integer, parameter        :: NWEATHER = 7
real                      :: MATRIX_WEATHER(NMAXDAYS,NWEATHER)
real   , dimension(100,3) :: CALENDAR_FERT, CALENDAR_NDEP, CALENDAR_PRUNT, CALENDAR_THINT
integer, dimension(100,2) :: DAYS_FERT    , DAYS_NDEP    , DAYS_PRUNT    , DAYS_THINT
real   , dimension(100)   :: NFERTV       , NDEPV        , FRPRUNT       , FRTHINT
real                      :: y(NDAYS,NOUT)

integer :: day, doy, NDAYS, NOUT, year

! State variables
real    :: chillday, Tsum
real    :: CR, CS, CB, CL, CRES, NL
real    :: WA, CLITT, CSOMF, CSOMS, NLITT, NSOMF, NSOMS, NMIN

! Additional output variables
real    :: ET_mmd, NPP_gCm2d, GPP_gCm2d, Reco_gCm2d, NEE_gCm2d


! PARAMETERS
call set_params(PARAMS)
! Composite parameters
NLAMAX   = NCLMAX / SLA     ! kg N m-2 leaf
NLAMIN   = NLAMAX * FNCLMIN ! kg N m-2 leaf

! ERRORS
call set_process_errors(PROCERR)
call set_state_errors(STATEERS)

! CALENDARS
! Calendar of weather
YEARI  = INT(MATRIX_WEATHER(:,1))
DOYI   = INT(MATRIX_WEATHER(:,2))
GRI    = MATRIX_WEATHER(:,3)
TI     = MATRIX_WEATHER(:,4)
RAINI  = MATRIX_WEATHER(:,5)
WNI    = MATRIX_WEATHER(:,6)
VPI    = MATRIX_WEATHER(:,7)
! Calendar of management
DAYS_FERT  = INT(CALENDAR_FERT (:,1:2))
DAYS_NDEP  = INT(CALENDAR_NDEP (:,1:2))
DAYS_PRUNT = INT(CALENDAR_PRUNT(:,1:2))
DAYS_THINT = INT(CALENDAR_THINT(:,1:2))
NFERTV     = CALENDAR_FERT (:,3)
NDEPV      = CALENDAR_NDEP (:,3)
FRPRUNT    = CALENDAR_PRUNT(:,3)
FRTHINT    = CALENDAR_THINT(:,3)


! Initial states
chillday = 0
Thist    = 10
! Thist    = TI(1) ! USE THIS CODE AFTER TESTING THE MODEL!!
Tsum     = 0

! Defining the inital parameters ( depending if model is restarted)

if (RESTART == 0) then
    treedens = TREEDENS0
    CR       = CRtree0 * TREEDENS0
    CS       = CStree0 * TREEDENS0
    CB       = CBtree0 * TREEDENS0
    CL       = CLtree0 * TREEDENS0
    NL       = CLtree0 * TREEDENS0 * NCLMAX
    if (FORTYPE == 2) then ! DECIDUOUS trees
      CRES   = CLtree0 * TREEDENS0 * 0.1
    else
      CRES   = 0.
    end if
    WA       = 1000 * ROOTD * WCST * FWCFC
    CLITT    = CLITT0
    CSOMF    = CSOM0 * FCSOMF0
    CSOMS    = CSOM0 * (1-FCSOMF0)
    NLITT    = CLITT0 / CNLITT0
    NSOMF    = (CSOM0 *    FCSOMF0)  / CNSOMF0
    NSOMS    = (CSOM0 * (1-FCSOMF0)) / CNSOMS0
    NMIN     = NMIN0

else
    treedens = STATEVARS(1)
    CR       = STATEVARS(2)
    CS       = STATEVARS(3)
    CB       = STATEVARS(4)
    CL       = STATEVARS(5)
    NL       = STATEVARS(6)
    CLITT = STATEVARS(7)
    CSOMF = STATEVARS(8)
    CSOMS = STATEVARS(9)
    NLITT = STATEVARS(10)
    NSOMF = STATEVARS(11)
    NSOMS = STATEVARS(12)
    NMIN = STATEVARS(13)
    WA = STATEVARS(14)

end if
do day = 1, NDAYS


! Environment
  call set_weather_day(day, year,doy)


! Tree phenology
  if ( FORTYPE == 2) then ! DECIDUOUS trees
    call dtsum_dchillday(FORTYPE,chillday,doy,Tsum, dchillday,dTsum)
    chillday = chillday + dchillday
    Tsum     = Tsum     + dTsum
  end if
  call phenology(FORTYPE,chillday,doy,Tsum, leaffall,treegrow)


! Management
  call fert_prune_thin(year,doy,DAYS_FERT,NFERTV,DAYS_PRUNT,FRPRUNT,DAYS_THINT,FRTHINT)
  treedens = treedens - thintreedens
  call morphology(CL,CS,CB,treedens,LAI)
  call foliarDynamics(CL,CRES,fTran,NL,LAI)


! Update model fluxes
  call PETtr(LAI,RAINint)
  call water_flux(WA,Evap,Tran,fTran,RWA,WFPS)
  call Nsupply(CR,NMIN,Nsup)
  call NPP(fTran, BIAS, RANDERR, fLUETERR)
  call allocation(FORTYPE,fTran,LAI, BIAS, RANDERR, FRERR)
  call NdemandOrgans
  call gtreeNupt(Nsup)
  call CNtree(CR,CS,CB)
  call water(WA,RAINint,Evap,Tran,LAI,RANDERR,BIAS,RUNOFFERR)
  call Tsoil_calc
  call CNsoil(RWA,WFPS,WA,gCR,CLITT,CSOMF,NLITT,NSOMF,NSOMS,NMIN,CSOMS, &
	BIAS, RANDERR, RESPERR, RLEACHERR)
  call N_dep(year,doy,DAYS_NDEP,NDEPV)



! Update model states
  CR    = CR    + gCR - dCR
  CS    = CS    + gCS - dCS
  CB    = CB    + gCB - dCB
  if (CL < 1.E-100) then
    CL  = 0.
    NL  = 0.
  else
    CL  = CL    + gCL - dCL + dCRESgrow
    NL  = NL    + gNL - dNL - retrNL + dCRESgrow * NCLMAX
  end if
  CRES = CRES + gCRES - dCRESgrow - dCRESthin
  WA    = WA    + RAIN - RAINint - Runoff - Drain - Evap - Tran
  CLITT = CLITT + dCL + dCB - rCLITT - dCLITT
  CSOMF = CSOMF + dCLITTsomf + dCR - rCSOMF - dCSOMF
  CSOMS = CSOMS + dCSOMFsoms - dCSOMS
  NLITT = NLITT + dNLlitt + dNBlitt - rNLITT - dNLITT
  NSOMF = NSOMF + dNRsomf + NLITTsomf - rNSOMF - dNSOMF
  NSOMS = NSOMS + NSOMFsoms - dNSOMS
  NMIN  = (NMIN  + Ndep + Nfert + Nmineralisation + Nfixation &
	- Nupt - Nleaching - Nemission)
  NMIN  = max(0.,NMIN)
	if (STATESPACE == 1) then
		NSOMS = NSOMS * NSOMSERR
		WA = WA*WAERR
		NL = NL*NPPERR
	end if
! Outputs
  y(day, 1)  = year + (doy-0.5)/366                   ! "Time" = Decimal year (approximation)
  y(day, 2)  = year
  y(day, 3)  = doy
  y(day, 4)  = treedens
  y(day, 5)  = CR
  y(day, 6)  = CS
  y(day, 7)  = CB
  y(day, 8)  = CL
  y(day, 9)  = NL
  y(day,10)  = CLITT
  y(day,11)  = CSOMF
  y(day,12)  = CSOMS
  y(day,13)  = NLITT
  y(day,14)  = NSOMF
	y(day,15)  = NSOMS
  y(day,16)  = NMIN
  y(day,17)  = WA
  NPP_gCm2d  = (gCL + gCB + gCS + gCR + gCRES) * 1000 ! g C m-2 d-1
  GPP_gCm2d  = NPP_gCm2d / (1-GAMMA)                  ! g C m-2 d-1
  Reco_gCm2d = (GPP_gCm2d - NPP_gCm2d) + Rsoil*1000   ! g C m-2 d-1
  NEE_gCm2d  = Reco_gCm2d - GPP_gCm2d                 ! g C m-2 d-1
  ET_mmd     = Evap + Tran                            ! mm d-1
  y(day,18)  = NEE_gCm2d                              ! g C m-2 d-1
  y(day,19)  = GPP_gCm2d                              ! g C m-2 d-1
  y(day,20)  = Reco_gCm2d                             ! g C m-2 d-1
  y(day,21)  = ET_mmd                                 ! mm d-1
  y(day,22) = WAERR
  y(day,23) = NSOMSERR
	y(day,24) = NPPERR
end do ! end time loop
end
