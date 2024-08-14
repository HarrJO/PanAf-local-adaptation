!17/01/2015: passage à la command line
!02/02/2015: introduction de la possibilité d'avoir des missing data dans le cas de données Poolseq seulement
!24/07/2015: passage à bypass avec modifs marginales (i) BF en dB, ii) resolution des Inf dans BPis, iii) retrait impression summary_omega si -omegafile, iv) ajout delta0yij option, v) si aux on active automatique covmcmc,)
!8/12/2015: corrections fuites memoires (cur_stream=0 qd pas dans boucle parallele, initialisation de dum_char_array la premiere fois), correction lecture opt_reg, opt_pheno et opt_aux avec initialisation propre, modification affichage ecran (anglais, matrice plus affichée, initialisation stream plus affichee, heure..) 

include 'mcmc_mat_utils.f90'
include 'updates_baypass.f90'
include 'M_kracken.f90' ! ligne de commande

program baypass
use mt_stream
use omp_lib
use mcmc_utils !fonctions diverses pour les summary statistics
use updates_baypass !update
use M_kracken

implicit none

integer, parameter :: max_line_length=100000 
! determine the maximal line length when reading geno input file (to count number of individual or population): on met un warning si plus de 10000 individuas car dans gl,
! si 10 caracteres par indiv. (incluant virgule et espace et phred like sur 2 chiffres)=> 10000 individus max
! augmenter la valeur=>augmentation du temps de parsing
integer:: err, i_thin,npop, nmrk, pop, mrk, nind, ind,max_nind, tst, tst_pi,tst_pij,tst_b,tst_cj,iter, pilot,&
          seed, nvaleurs, thin,burn_in,npilot,pilot_length,y_out,tst_y,&
          dum_int,dum_int2,ddl_wish,rho,tmp_y_rep,npheno,pheno,ninter_beta,delta0_yij,&
          iflen,igot,ibegin(2),iterm(2),ilen,ier,nstream,cur_stream,root_poolsize,contrast,ncontrast
integer (kind=8) ::  val_systime(3)         
integer, allocatable :: Y_OBS(:,:), N_OBS(:,:) ,& !observations
                        PoolSize(:),Y_READS(:,:),N_READS(:,:),delta_y(:,:),Yminmax(:,:,:),&
                        INITS_Delta(:,:),&
                        R_CNT(:,:),R_READS(:,:),delta_rooty(:),RootYminmax(:,:),& !R_CNT(nmrk,3): pour la racine, cnt_all1,cnt_all2,cnt_tot: initialise à zero par defaut car pas d'info sur la racine; si root=pool => R_READS et R_Delta
                        contrastes(:,:),& !ncontrast,npop: doit etre -1,0 ou 1
                        omega_index_subpop(:,:),dum_int_pop(:),& !tab de dim npop,npop-1 dans lesquels sont stockés les indices relevant pour les extractions de chque pop (cf alpha_up_muvar)
                        dum_read_line(:),ind_pop_index(:),nind_par_pop(:)
real (kind=8), allocatable :: INITS_Pij(:,:,:),mean_pij(:,:,:),delta_pij(:,:),acc_pij(:,:),& !frequences alleliques
                              INITS_Pi(:),delta_pi(:),acc_pi(:),mean_pi(:,:),& !frequences ancestrales
                              delta_nichc(:),acc_nichc(:),&! si option nicholson
                              PHENO_VAL(:,:),INITS_Beta_i(:,:),INITS_tau_beta(:),INIT_Pdelta(:),&
                              mean_beta_i(:,:,:),mean_tau_beta(:,:),mean_delta_i(:,:),mean_p_delta(:,:),sum_covariates(:,:),&
                              delta_betai(:,:),acc_betai(:,:),&!au cas ou up_betai_coop
                              INITS_LDA(:,:),mean_lda(:,:,:),omega_mat(:,:),mean_omega(:,:,:),mean_xtx(:,:),mean_contrast(:,:,:),contrast_sd(:),& 
                              acc_y(:,:),acc_rooty(:),CPO(:,:),log_cn(:,:),mean_yij(:,:,:),mean_rooty(:,:),muvar_cur(:,:),&
                              dum_real_pop(:),dum_real_pop2(:),dum_real_pheno(:),c_mat(:,:),dum_vect_grid(:,:),&
                              reg_mat(:,:),inv_c_mat(:,:),mean_rho_coef(:,:,:),mean_isbf(:,:,:),grid_is(:,:,:),dum_pheno_std(:,:),dum_freq(:,:,:),&
                              array_real_pop(:,:,:),PHENO_OBS(:,:,:),mean_phenoval(:,:,:),& !pour pheno_obs: 1=sd à la lecture puis transformation en precision (0 si missing)
                              cur_freq_like(:,:),ind_pl(:,:,:,:) !genotype likelihood de dimensions (nmrk,npop,nmaxind_par_pop,3); 
                              !Allele frequency are with respect to the reference allele); the four last dimension elements are i) PL homo ref (0); ii) PL hetero (1); iii) PL homo alt (2)                              
real (kind=8):: delta0_pij,delta0_pi,delta0_cj,acc_inf=0.25,acc_sup=0.4,rate_adjust=1.25,niter_tot,& !pour etudes pilotes
             !  beta_pi=0.7,lda_a_prior=3.,lda_b_prior=4.,& !parmatres hyperpriors (a rentrer)
                pi_beta_params(2),acc_beta_params(2),mean_beta_params(2,2),pi_beta_update_save(5),& !Sum(log(Pi)) ; Sum(log(1-pi)) ; gamma_log(Phi) ; gamma_log(Phi*mu) ; gamma_log(Phi*(1-mu)) avec Phi=(a_pi+b_pi) et mu=a_pi/Phi
                P_beta_params(2),&
                delta_pi_beta_mu=0.05,delta_pi_beta_phi=1,&
                bpval,kld,bf,tmp_mean,tmp_min,tmp_max,dum_real,dum_real2,dum_real5(5),tmp_a,tmp_b,&
                alpha_up_muvar_param(2),tmp_t_rep,tmp_t,tmp_t_obs,tmp_vij,k_tau_prior=1.,l_tau_prior=10.,tau_beta0=20.,&
                min_beta=-0.3 , max_beta=0.3 ,ising_b,covgaussprior(3) !covgaussprior(mean,tau,sd)
real (kind=8) :: dum_real16,dum_real16_2, mean_deviance,lpml !passage a kind=8 car pb avc certains compulateur
!attention tres important d'initialiser ici les opt et dum_char_array
logical  :: up_params_beta, poolseq=.false. , up_alpha_coop=.true., opt_pheno=.false. ,& !up_slc_pi=.false.
            opt_aux=.false. , estim_tau_beta=.false. , out_pilot , estim_omega,opt_reg=.false.,&
            up_betai_coop,opt_ising,scale_cov,nich_prior=.false.,printintmat=.false.,printsampledomega=.false.,root_geno=.false.,root_pool=.false.,opt_contrast=.false.,&
            opt_predcov=.false.,indseq=.false.
logical, allocatable :: missing_data(:,:),root_missing(:)

character(len=255) :: Geno_file,Pheno_file,Pheno_SDfile,Root_Geno_file,omega_file,poolsize_file,contrast_file,dum_char='',dum_char_array(2)='',&
                      sum_yij_pij_file,sum_pij_file,sum_pi_file,sum_omega_file,sum_mat_omega_file,&
                      sum_betaparams_file,sum_betai_file,sum_Pdelta_file,sum_taufile,&
                      sum_beta_i_reg_file,sum_DIC_file,sum_phenostd_file,out_prefix,sum_covpred_file,sum_contrast_file,sampled_omega_file
character(len=max_line_length)  :: big_string                    

type(mt_state),allocatable :: mts(:) !mts(0:NSTREAM-1)

call system_clock(val_systime(1),val_systime(2),val_systime(3))
call void_read_input()
niter_tot=real(npilot)*real(pilot_length)+real(burn_in)+real(nvaleurs)*real(thin)
!~ call sgrnd(seed) !seed du mersenne twister

print *,''
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'        PILOT RUNS START'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,''

pilot=1 
 call void_initialisation()
 print *,'' ; print *,'###########INITS VALUES#############' ; print *,''
 call void_printval()

do while (pilot<=npilot)
 acc_pij=0. ; acc_pi=0.  ; acc_beta_params=0.
 if(nich_prior) acc_nichc=0.
 if(opt_pheno .and. up_betai_coop) acc_betai=0.
 if(poolseq) acc_y=0.
 if(root_pool) acc_rooty=0.
 print *,'PILOT RUN: ',pilot
 do iter=1,pilot_length
  if(mod(iter,pilot_length/10)==0) print *,'  iteration=',iter
  call void_mcmc_iter()
 end do
 call void_adjust_proposals()
 print *,'' ; print *,'########### NODE VALUES #############' ; print *,''
 call void_printval()
 pilot=pilot+1
 call void_initialisation()
 print *,''
 call elapse_time(val_systime,real(pilot)*real(pilot_length),niter_tot)
 print *,''
end do
call elapse_time(val_systime,-1.,-1.)

print *,''
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'        BURN-IN: START'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,''

if(opt_aux) then
 if(estim_tau_beta) then
  allocate(dum_real_pheno(npheno)) ;  dum_real_pheno=0.
 else
  out_pilot=.true.
 end if
end if

do iter=1,burn_in
 if(mod(iter,burn_in/10)==0) then
  print *,'  iteration=',iter
  call elapse_time(val_systime,real(pilot)*real(pilot_length)+real(iter),niter_tot)
 end if
 call void_mcmc_iter()
 if(opt_aux .and. estim_tau_beta) dum_real_pheno=dum_real_pheno + INITS_tau_beta
end do

if(opt_aux .and. estim_tau_beta) then
 INITS_tau_beta=(dum_real_pheno/burn_in)*(P_beta_params(1)/sum(P_beta_params))
 out_pilot=.true.
 print *,'' ; print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' ; print *,'     2nd BURN-IN: START'
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' ; print *,''
 write(*,'(A,1x,100(f12.6,1x))') 'Estimated Precision(s): ',INITS_tau_beta
 print *,''
 do iter=1,burn_in
  if(mod(iter,burn_in/10)==0) then 
   print *,'  iteration=',iter
   call elapse_time(val_systime,real(pilot)*real(pilot_length)+real(iter),niter_tot) ! a revoir
  end if
  call void_mcmc_iter()
 end do
end if

call void_printval()
call elapse_time(val_systime,-1.,-1.)
print *,''
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'        MCMC CHAIN STARTED'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,''
print *,''

 allocate(CPO(nmrk,npop)) ; CPO=0. ; mean_deviance=0.
 acc_pij=0. ; acc_pi=0. ; acc_beta_params=0.
 allocate(mean_pij(nmrk,npop,4),mean_pi(nmrk,2),mean_xtx(nmrk,4),mean_lda(npop,npop,2),mean_omega(npop,npop,2))
 mean_lda=0. ; mean_omega=0. ; mean_pij=0.; mean_pi=0.; mean_beta_params=0. ; mean_xtx=0.
 if(poolseq) then
  allocate(mean_yij(nmrk,npop,2))
  mean_yij=0. ; acc_y=0.
 end if

 if(root_pool) then
  allocate(mean_rooty(nmrk,2))
  mean_rooty=0. ; acc_rooty=0.
 end if

 if(nich_prior) acc_nichc=0.

 if(opt_pheno) then
  allocate(mean_beta_i(nmrk,npheno,2),mean_tau_beta(npheno,2),mean_delta_i(nmrk,npheno),mean_p_delta(npheno,2))
  mean_beta_i=0. ; mean_tau_beta=0. ; mean_delta_i=0. ; mean_p_delta=0.
  if(opt_predcov) then
   allocate(mean_phenoval(npop,npheno,2))
   mean_phenoval=0.
  end if
 end if

 if(printsampledomega) then
  open(101,file=sampled_omega_file,status='unknown')
!  write(101,*) 'Npops: ',npop
 end if 
 do iter=1,nvaleurs
  if(nvaleurs>100 .and. mod(iter,nvaleurs/100)==0) then
   print *,'  iteration=',iter
   call elapse_time(val_systime,real(pilot)*real(pilot_length)+real(burn_in)+real(iter)*real(thin),niter_tot)
  end if
  do i_thin=1,thin
   call void_mcmc_iter()
  end do
  if(printsampledomega) then
  do pop=1,npop
   write(101,'(1000(f12.8,1x))') ( omega_mat(pop,dum_int), dum_int=1,npop )
  end do
  end if
  call void_update_summary()
 end do
 if(printsampledomega) close(101)
 call void_print_summary()

 do pop=0,nstream-1
    call delete(mts(pop))
 enddo 

print *,''
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'        END '
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,''
call elapse_time(val_systime,-1.,-1.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! calcul des differentes reg_mat (il y en a npop)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_reg_mat() 
 do dum_int=1,npop
  reg_mat(dum_int,:)=matmul(omega_mat(dum_int,omega_index_subpop(dum_int,:)),inv(omega_mat(omega_index_subpop(dum_int,:),omega_index_subpop(dum_int,:))))
 end do
end subroutine compute_reg_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   UPDATE DE MATRICE LDA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function up_lda_mat()
 integer (kind=8) :: tmp_i,tmp_j
 real (kind=8) :: up_lda_mat(npop,npop),mat_wish(npop,npop) 

  mat_wish=0.
  if(opt_pheno) then
   do tmp_i=1,nmrk
    mat_wish=mat_wish + txx(INITS_Pij(tmp_i,:,2) - sum_covariates(tmp_i,:)/(sqrt(INITS_Pi(tmp_i)*(1.-INITS_Pi(tmp_i)))))
   end do
  else
   do tmp_i=1,nmrk
    mat_wish=mat_wish + txx(INITS_Pij(tmp_i,:,2))
   end do
  end if
  do tmp_j=1,npop
    mat_wish(tmp_j,tmp_j)=mat_wish(tmp_j,tmp_j) + rho
  end do
  mat_wish=inv(mat_wish)
  call wishart_sample ( npop, ddl_wish, mat_wish, up_lda_mat, mts(cur_stream) )

 end function up_lda_mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!LECTURE DES INPUT/ALLOCATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_read_input()
 call kracken('cmd','-accinf 0.25 -accsup 0.4 -adjrate 1.25 -auxmodel no -auxPbetaprior 0.02 1.98 -betapiprior 1.0 1.0 -burnin 5000 -contrastfile -covmcmc no -covpriormean 0. -covpriorsd 1. -d0cj 0.05 -d0pi 0.5 -d0pij 0.05 -d0yij 1 -efile -efilesd -setpibetapar no -esttaubeta no -gfile -help no -indseq no -isingbeta 0.0 -maxbeta 0.3 -minbeta -0.3 -nbetagrid 201 -nicholsonprior no -npilot 20 -npop -1 -nthreads 1 -nval 1000 -omegafile -outprefix -pilotlength 500 -printintmat no -print_omega_samples no -poolsizefile -rho 1 -rootgfile -rootpoolsize -1 -scalecov no -seed 5001 -tauB0 -1. -thin 20 -upalphaalt no')      

 call retrev("cmd_help",dum_char,iflen,ier)
 if(iflen==0) then !on imprime l'aide
  call void_print_help()
  stop
 end if

 call retrev("cmd_indseq",dum_char,iflen,ier)
 if(iflen==0) indseq=.true.
!~  npop = iget("cmd_npop")
!~  if(npop<0) then
!~   write(*,'(A)') 'ERROR: Please provide the number of populations'
!~   stop
!~  end if

!recuperation de nmrk et npop et on verifie que ca s'ouvre bien
 call retrev("cmd_gfile",dum_char,iflen,ier)
 if(iflen==0) then
  write(*,'(A)') 'ERROR: Please provide a genotyping file'
 stop
 end if
 Geno_file=dum_char(:iflen)
 open(1,file=Geno_file,status='old') 
 if(indseq) then !il y a alors un en-tete qui contient les index des pops pour chaque individu
  read(1,'(a)',iostat=err) big_string
  allocate(dum_read_line(max_line_length))
   do ind=1,len_trim(big_string) 
    read(big_string, *, iostat=err) dum_read_line(1:ind)
   if(err/=0) exit
   end do
   nind=ind-1
   if(nind>10000) then
    write(*,'(A)') 'WARNING: More than 10000 individuals: this may cause failure in parsing the genofile'
    write(*,'(A)') 'WARNING: Either subsample the number of ind or recompile after increasing max_line_length value in the baypass.f90 source file (L19)'    
   end if
   allocate(ind_pop_index(nind))
   ind_pop_index=dum_read_line(1:nind) 
   deallocate(dum_read_line)
   !count number of individual per pop and check that all pops are represented by at least one individual
   npop=maxval(ind_pop_index)
   if(npop<2) then
    write(*,'(A)') 'ERROR: At least 2 populations must be represented'
   end if
   allocate(nind_par_pop(npop))
   nind_par_pop=0
   do ind=1,nind
     nind_par_pop(ind_pop_index(ind))=nind_par_pop(ind_pop_index(ind))+1
   end do
   !check nombre d'individu par pop
    do pop=1,npop
     if(nind_par_pop(pop)==0) then
       write(*,'(A)') 'ERROR: No individual for pop Num. ',pop
       write(*,'(A)') 'Check Individual assignation (first line of genopl file)'
       stop
     end if
    end do
   up_alpha_coop=.false. !on est dans le meme cas que update_alpha_alternatif (i.e., freq par freq plutot que vecteur) => permet d'avoir les bonnes dimensions pour les arrays des parametres mcmc (acc, delta) et de faire les decompositions de omega 
 end if  
 
 !decompte des nombre de marqueurs
 nmrk=0 !Rk si indseq l'en tete vient d etre lue au dessus
 do 
  read(1,*,end=1,iostat=err) 
  nmrk=nmrk+1
 end do
 1 continue
 close(1)
 
 !decompte du nombre de pop sur le fichier geno et verification toutes les lignes ont le meme nombre de pops
 if(.not. indseq) then !count number of pop
  open(1,file=Geno_file,status='old') 
  read(1,'(a)',iostat=err) big_string
  allocate(dum_read_line(max_line_length))
   do pop=1,len_trim(big_string) 
    read(big_string, *, iostat=err) dum_read_line(1:pop)
   if(err/=0) exit
   end do
   npop=(pop-1)/2 
   if(npop<2) then
    write(*,'(A)') 'ERROR: At least 2 populations must be represented'
    stop
   end if
 !verif sur le reste du fichier qu'on a toujours le meme nombre de pop estimees
  allocate(dum_real_pop(2*npop+1))
  do mrk=1,nmrk-1
  read(1,'(a)',iostat=err) big_string
  read(big_string, *, iostat=err) dum_real_pop(1:2*npop) 
  read(big_string, *, iostat=dum_int) dum_real_pop 
   if(err/=0 .or. dum_int==0) then
    write(*,'(A,i8,A)') 'ERROR in genotype file (declared with -gfile) for line number',mrk+1,' :'
    write(*,'(A)') trim(big_string)
    write(*,'(A,i4,A,i4,A)') 'It should contain ',2*npop,' elements (since ',npop,' pops were inferred from the number of elements in the first line but line)'
    stop
   end if 
  end do
  deallocate(dum_read_line,dum_real_pop)  
  close(1)  
 end if

!recuperation eventuelle info sur la pop racine (et test que le fichier geno a le meme nombre de ligne
 allocate(R_CNT(nmrk,3)) ! tres important d'initialiser car utilise pour update des pi (doit etre a 0 si pas d'info sur la pop root)
 R_CNT=0 
 call retrev("cmd_rootgfile",dum_char,iflen,ier) 
 if(iflen/=0) then
  root_geno=.true.
  Root_Geno_file=dum_char(:iflen)
  allocate(dum_real_pop(3)) !check du nombre d'elements: que 2
  open(1,file=Root_Geno_file,status='old') 
  dum_int=0
  do 
   read(1,'(a)',end=2,iostat=err) big_string
   dum_int=dum_int+1
   read(big_string, *, iostat=err) dum_real_pop(1:2) 
   read(big_string, *, iostat=dum_int2) dum_real_pop 
   if(err/=0 .or. dum_int2==0) then
    write(*,'(A,i8,A)') 'ERROR in root genotype file (declared with -rootgfile) for line number',dum_int,' :'
    write(*,'(A)') trim(big_string)
    write(*,'(A)') 'It should contain two elements (counts for all1 and counts for allele 2 in the root population)'
    stop
   end if       
  end do
  2 continue
  deallocate(dum_real_pop)
  close(1)
  if(dum_int/=nmrk) then
   write(*,'(A)') 'Genotyping files for the root population and the current populations have different numbers of markers'
   stop  
  end if
  root_poolsize=iget("cmd_rootpoolsize") !sert de test pour savoir si racine=pool (si<0=>comptages sinon reads
  if(root_poolsize>0) root_pool=.true.
 end if


 seed=iget("cmd_seed") ; ninter_beta=iget("cmd_nbetagrid")
 npilot = iget("cmd_npilot") ; burn_in =  iget("cmd_burnin") ; thin=iget("cmd_thin")
 nvaleurs=iget("cmd_nval") ; pilot_length = iget("cmd_pilotlength")
 delta0_cj =rget("cmd_d0cj") ; delta0_pi=rget("cmd_d0pi") ; delta0_pij=rget("cmd_d0pij") ; delta0_yij=iget("cmd_d0yij")
 rate_adjust=rget("cmd_adjrate") ; acc_inf=rget("cmd_accinf") ; acc_sup=rget("cmd_accsup")
 min_beta=rget("cmd_minbeta") ; max_beta=rget("cmd_maxbeta")

 write(*,'(A,i9)') ' No of Markers                      = ',nmrk 
 if(indseq) then
  write(*,'(A)') ' INDSEQ mode:   '  
  write(*,'(A,i9)') ' No of Individuals                  = ',nind
  write(*,'(A,i9)') ' No of Populations                  = ',npop   
  tmp_mean=0.;dum_int=nind;max_nind=0
  do pop=1,npop
   tmp_mean=tmp_mean+nind_par_pop(pop)
   if(nind_par_pop(pop)>max_nind) max_nind=nind_par_pop(pop)
   if(nind_par_pop(pop)<dum_int)  dum_int=nind_par_pop(pop) 
 end do
 write(*,'(t3,A,i5,A,i5,A,f6.2,A)') ' From ', dum_int,' to ',max_nind,' ind. per pop (mean=',tmp_mean/npop,')'
 else
  write(*,'(A,i9)') ' No of Populations                  = ',npop
 end if
 write(*,'(A,A)') ' Genotype File name                 =    ',trim(Geno_file)
 print *,''
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,'        Model Specifications'
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,''

!poolseq?
 call retrev("cmd_poolsizefile",dum_char,iflen,ier) 
 if(iflen>0) then
  if(indseq) then
   write(*,'(A)') 'ERROR: Option incompatibilities: indseq mode is not compatible with PoolSeq mode: Check options'
  end if
  poolseq=.true. ; poolsize_file=dum_char(:iflen)
  allocate(PoolSize(npop),dum_real_pop(npop+1))
  poolsize_file=dum_char(:iflen)
  open(1,file=poolsize_file,status='old')
  read(1,'(a)',iostat=err) big_string
  read(big_string, *, iostat=err) PoolSize
  read(big_string, *, iostat=dum_int) dum_real_pop 
   if(err/=0 .or. dum_int==0) then
    write(*,'(A)') 'ERROR: The number of elements in the poolsize file (declared with -poolsizefile) is not correct:'
    write(*,'(A)') trim(big_string)
    write(*,'(A,i4,A,i4,A)') 'It should contain ',npop,' elements (since ',npop,' pops were inferred from the geno file)'
    stop
   else
    write(*,'(A,1000(i5,1x))') ' -> POOLSEQ MODE with pool sizes: ',PoolSize
   end if   
  close(1)
  deallocate(dum_real_pop) 
 end if

!beta prior sur les Pi?
 call retrev("cmd_setpibetapar",dum_char,iflen,ier)
 up_params_beta=.true. ; if(iflen==0) up_params_beta=.false.
 print *,''
 call retrev("cmd_betapiprior",dum_char,iflen,ier)
 call delim(dum_char,dum_char_array(1:2),2,igot,ibegin,iterm,ilen,' ,:')
 call string_to_real(dum_char_array(1),pi_beta_params(1),ier)
 call string_to_real(dum_char_array(2),pi_beta_params(2),ier)

!!Type de modele possible: i)   sans covariate (avec ou sans estim_omega)  =>opt_pheno=opt_reg=opt_aux=false et estim_omega=true ou false
!                          ii)  IS covariate (avec ou sans estim_omega)     =>opt_reg=true ; opt_pheno=opt_aux=false   et estim_omega=true ou false
!                          iii) MCMC covariate (forcement sans estim_omega) => STD model (opt_pheno=true et opt_aux=false) OU AUX model (opt_pheno=true et opt_aux=true);  et estim_omega=false,opt_ppp=.false
!                               iii.1) prior uniforme sur beta avec ou sans aux
!                               iii.2) prior gaussienne avec ou sans aux

!on recupere l'ensemble des infos et on gere les incompatibilites'

 call retrev("cmd_omegafile",dum_char,iflen,ier)  
 if(iflen>0) then
  omega_file=dum_char(:iflen)
  estim_omega=.false.
 else 
  estim_omega=.true.
  call retrev("cmd_nicholsonprior",dum_char,iflen,ier)
  if(iflen==0) then
    nich_prior=.true.
    allocate(delta_nichc(npop),acc_nichc(npop))
    delta_nichc=delta0_cj
  else
   rho=iget("cmd_rho") ; ddl_wish=iget("cmd_rho") + nmrk 
  end if
 end if
 call retrev("cmd_printintmat",dum_char,iflen,ier)
 if(iflen==0) printintmat=.true. !on imprime omega au cours des etudes pilotes la matrice
 call retrev("cmd_print_omega_samples",dum_char,iflen,ier)
 if(iflen==0) printsampledomega=.true. !on imprime omega à chaque post-burnin et thinned MCMC sample
 
 call retrev("cmd_auxPbetaprior",dum_char,iflen,ier)
 call delim(dum_char,dum_char_array(1:2),2,igot,ibegin,iterm,ilen,' ,:')
 call string_to_real(dum_char_array(1),P_beta_params(1),ier)
 call string_to_real(dum_char_array(2),P_beta_params(2),ier)

 call retrev("cmd_esttaubeta",dum_char,iflen,ier)
 estim_tau_beta=.false. ; if(iflen==0) estim_tau_beta=.true.
 tau_beta0=rget("cmd_tauB0") 
 if(tau_beta0<0) then
  up_betai_coop=.true.
 else
  up_betai_coop=.false.
 end if

!!!!!!!!!!!!!!!!!!!!
!!!les modeles:
!!!!!!!!!!!!!!!!!!!!

 call retrev("cmd_efile",dum_char,iflen,ier) 
 !ouvrir et recuperer npheno (on verifie par la meme occasion que le nombre d'elements est OK dans chaque ligne)
 if(iflen>0) then
  opt_reg=.true. ! par defaut modele IS
  Pheno_file=dum_char(:iflen) !recuperation de npheno
  allocate(dum_real_pop(npop+1)) !pour controler le nombre d'elements
  open(1,file=Pheno_file,status='old') 
  npheno=0
  do 
   read(1,'(a)',end=3,iostat=err) big_string
   npheno=npheno+1
   read(big_string, *, iostat=err) dum_real_pop(1:npop)
   read(big_string, *, iostat=dum_int) dum_real_pop 
   if(err/=0 .or. dum_int==0) then
    write(*,'(A,i3)') 'ERROR: The number of elements in the covariate file (declared with -efile) is not correct for line: ',npheno
    write(*,'(A)') trim(big_string)
    write(*,'(A,i4,A,i4,A)') 'It should contain ',npop,' elements (since ',npop,' pops were inferred from the geno file)'
    stop   
   end if
  end do
  3 continue
  deallocate(dum_real_pop) 
  close(1)
  call retrev("cmd_scalecov",dum_char,iflen,ier)
  scale_cov=.false. ; if(iflen==0) scale_cov=.true.
 end if

!is ou mcmc mode?
 call retrev("cmd_covmcmc",dum_char,iflen,ier)
 if(iflen==0 .or. opt_aux) then
  opt_pheno=.true.
  if(opt_reg) then !on avait bien un fichier efile mais on veut desactiver IS
   opt_reg=.false.
  else
    print *,'The covmcmc mode requires a covariate file (-efile option)'
    stop
  end if
 end if

 call retrev("cmd_auxmodel",dum_char,iflen,ier)
 if(iflen==0) then
  opt_aux=.true. ;  opt_pheno=.true.
  if(opt_reg) then !on avait bien un fichier efile mais on veut desactiver IS
   opt_reg=.false.
  else
    print *,'The covmcmc mode requires a covariate file (-efile option)'
    stop
  end if
 end if
 ising_b=rget("cmd_isingbeta") 

 if(opt_pheno) then
  if(estim_omega) then
    print *,'ERROR: MCMC STD or AUX covariate mode is activated but no omega matrix file was provided'
    print *,' Please provide an Omega matrix'
    print *,' AND/OR switch to the IS covariate mode (i.e. remove -covmcmc or -aux option)'
    stop
   end if
 end if

!Covariate prediction mode
 call retrev("cmd_efilesd",dum_char,iflen,ier)
 if(iflen>0) then
  if(.not. opt_pheno) then
    print *,'ERROR: Prediction Covariate mode requires activation of the MCMC STD or AUX covariate mode'
    print *,'Please activate -covmcmc or -aux option (and provide an omegafile)'
    stop  
  end if
  opt_predcov=.true.
  Pheno_SDfile=dum_char(:iflen) 
  dum_int2=0 ! controle du nombre de ligne
  allocate(dum_real_pop(npop+1)) !pour controler le nombre d'elements
  open(1,file=Pheno_SDfile,status='old')
  do 
   read(1,'(a)',end=4,iostat=err) big_string
   dum_int2=dum_int2+1
   read(big_string, *, iostat=err) dum_real_pop(1:npop)
   read(big_string, *, iostat=dum_int) dum_real_pop 
   if(err/=0 .or. dum_int==0) then
    write(*,'(A,i3)') 'ERROR: The number of elements in the sd covariate file (declared with -efilesd) is not correct for line: ',dum_int2
    write(*,'(A)') trim(big_string)
    write(*,'(A,i4,A,i4,A)') 'It should contain ',npop,' elements (since ',npop,' pops were inferred from the geno file)'
    stop   
   end if
  end do
  4 continue
  deallocate(dum_real_pop)
  close(1)
  if(dum_int2/=npheno) then
   print *,'ERROR: Covariate SD file and Covariate have different raw (covariate) numbers'
   stop
  end if
  covgaussprior(1) = rget("cmd_covpriormean") ; covgaussprior(3) = rget("cmd_covpriorsd") ; covgaussprior(2) = 1./(covgaussprior(3)*covgaussprior(3))
 end if

!Calcul des contrastes (le cas echeant)
 call retrev("cmd_contrastfile",dum_char,iflen,ier) 
 !ouvrir et recuperer ncontrast 
 if(iflen>0) then
  opt_contrast=.true. 
  contrast_file=dum_char(:iflen) !recuperation du nombre de contraste et verification qu'il est bien egal au nombre de pops
  allocate(dum_real_pop(npop+1)) !pour controler le nombre d'elements
  open(1,file=contrast_file,status='old') 
  ncontrast=0
  do 
   read(1,'(a)',end=5,iostat=err)  big_string
   ncontrast=ncontrast+1
   read(big_string, *, iostat=err) dum_real_pop(1:npop)
   read(big_string, *, iostat=dum_int) dum_real_pop 
   if(err/=0 .or. dum_int==0) then
    write(*,'(A,i3)') 'ERROR: The number of elements in the contrast file (declared with -contrastfile) is not correct for line: ',ncontrast
    write(*,'(A)') trim(big_string)
    write(*,'(A,i4,A,i4,A)') 'It should contain ',npop,' elements (since ',npop,' pops were inferred from the geno file)'
    stop   
   end if
  end do
  5 continue
  deallocate(dum_real_pop)
  close(1)
 end if

!recuperation du prefixe output
 call retrev("cmd_outprefix",dum_char,iflen,ier)
 if(iflen==0) then
  out_prefix=''
 else 
   out_prefix=trim(dum_char(:iflen))//'_'
 end if 
 sum_yij_pij_file=trim(out_prefix)//'summary_yij_pij.out' ;  sum_pij_file=trim(out_prefix)//'summary_pij.out'
 sum_pi_file=trim(out_prefix)//'summary_pi_xtx.out'  ; sum_omega_file=trim(out_prefix)//'summary_lda_omega.out'
 sum_mat_omega_file=trim(out_prefix)//'mat_omega.out' ; sum_betaparams_file=trim(out_prefix)//'summary_beta_params.out'
 sum_betai_file=trim(out_prefix)//'summary_betai.out' ; sum_Pdelta_file=trim(out_prefix)//'summary_Pdelta.out'
 sum_taufile=trim(out_prefix)//'summary_tau.out' ; sum_beta_i_reg_file=trim(out_prefix)//'summary_betai_reg.out'
 sum_DIC_file=trim(out_prefix)//'DIC.out' ; sum_phenostd_file=trim(out_prefix)//'covariate.std'
 sum_covpred_file=trim(out_prefix)//'predcov.out' ; sum_contrast_file=trim(out_prefix)//'summary_contrast.out'
 sampled_omega_file=trim(out_prefix)//'omegasamples.out'

!recuperation des threads et initialisation du RNG
 nstream=iget("cmd_nthreads")
 call omp_set_num_threads ( nstream )
 call set_mt19937
 allocate(mts(0:(nstream-1)))
 call new(mts(0))
 call init(mts(0),seed)  ! init by scalar
 do pop=1,nstream-1
    call create_stream(mts(0),mts(pop),pop)
 enddo 

!!!!!!!!!!
!Synthese:
!!!!!!!!!!

!modele basique
 if(.not. opt_pheno .and. .not. opt_reg) then
  print *,' -> Basic Mode (no covariate):'
  if(estim_omega) then
   if(nich_prior) then
    write(*,'(A)') '  (*) Omega matrix is estimated assuming a Nicholson prior (i.e.: Omega is diagonal with cj~Unif(0,1)' 
   else
    write(*,'(A,i5)') '  (*) Omega matrix is estimated with prior inv(Omega) ~ Wish_npop((1/rho)Id_npop,rho)  where rho= ',rho
   end if
  end if
  if(up_params_beta) then 
   write(*,'(A)') '  (*) Parameters from the Pi prior dist. will be estimated'
  else
   write(*,'(A,f6.3,A,f6.3,A)') '  (*) Prior on Pi is Beta(',pi_beta_params(1),',',pi_beta_params(2),')'
  end if
   write(*,'(A)') '  (*) Outlier detection: XtX stat. for each SNP        '
 end if

!calcul des contrastes
 if(opt_contrast) then
  write(*,('(A,i5)')) '  (*) Contrasts between std. allele frequencies will be computed: ncontrast = ',ncontrast
 end if

!modele IS covariate
 if(opt_reg) then
   write(*,'(A,i3,A)') '  -> IS Covariate Mode with ',npheno, ' covariates (i.e. BF and Beta coef will be estimated using an Importance Sampling algorithm):'
  if(estim_omega) &
   write(*,'(A,i5)') '  (*) Omega matrix is estimated with prior inv(Omega) ~ Wish_npop((1/rho)Id_npop,rho)  where rho= ',rho
  if(up_params_beta) then 
   write(*,'(A)') '  (*) Parameters from the Pi prior dist. will be estimated'
  else
   write(*,'(A,f6.3,A,f6.3,A)') '  (*) Prior on Pi is Beta(',pi_beta_params(1),',',pi_beta_params(2),')'
  end if
   write(*,'(A)') '  (*) Outlier detection: XtX stat. for each SNP        '
   write(*,('(A,f8.4,A,f8.4,A)')) '  (*) Prior on Beta coef is U(',min_beta,',',max_beta,')'
   write(*,('(A,i4,A,f8.4,A,f8.4,A)')) '  (*) IS estimates are approximate on the grid with ',ninter_beta,' points uniformly distributed over the (',min_beta,',',max_beta,') interval' 
 end if

!modele MCMC covariate
 if(opt_pheno) then
  write(*,'(A,i3,A)') '  -> MCMC Covariate Mode with ',npheno, ' covariates (outlier detection and estimation of omega are inactivated)'
   if(up_betai_coop) then
    write(*,('(A,f8.4,A,f8.4,A)')) '  (*) Prior on Beta coef is U(',min_beta,',',max_beta,')'
   else
    write(*,('(A,f8.4,A,f8.4,A)')) '  (*) Prior on Beta coef is N(0,',1./tau_beta0,')'
   end if
   if(opt_aux) then
    write(*,('(A,f6.3,A,f6.3,A)')) '  (*) Auxiliary variable model with P~Beta(',P_beta_params(1),',',P_beta_params(2),') => BF will be calculated'
    if(ising_b>-0.0001 .and. ising_b<0.0001) then
     write(*,'(A)') '        => No spatial dependancy between markers (isingbeta=0.)'
     ising_b=0.
     opt_ising=.false.
    else 
     opt_ising=.true.
     write(*,('(A,f6.3)')) '        => Spatial dependancy between markers is modeled via an Ising model on the auxiliary variable with beta_ising=',ising_b
    end if
   else
    write(*,'(A)') '  (*) Full model (BPvalue will be computed)'
   end if
   if(opt_predcov) then
    write(*,('(A,f6.3,A,f6.3,A)')) '  (*) Covariate Prediction Mode assuming a prior Zj~N(',covgaussprior(1),',',covgaussprior(3),') => Missing pop. covariable values will be predicted'     
   end if
 end if

 print *,''
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,'           MCMC specifications'
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,''
 print *,' Nb. of sampled parameter values       = ',nvaleurs
 print *,' Thinning Rate                         = ',thin
 print *,' Burn in Period Length                 = ',burn_in
 print *,' Max Number of Pilot runs              = ',npilot
 print *,' Pilot run Length                      = ',pilot_length
 print *,''

 write(*,'(A,1x,f6.4)')  ' Init. deltas for Pi proposals          = ',delta0_pi
 call retrev("cmd_upalphaalt",dum_char,iflen,ier)
 if(iflen==0) up_alpha_coop=.false. !indseq controle ensuite si update comptage ou genotype likelihood
 if(up_alpha_coop) then
  write(*,'(A,1x,f6.4)') ' Coop Pij proposal with init. sig. for Pi. multivariate prop.          = ',delta0_pij 
 else
  write(*,'(A,1x,f6.4)') ' Pij proposal with init. sig. for Pij prop.                       = ',delta0_pij
 end if
 write(*,'(A,1x,2(f4.2,1x))') ' Pilot Run Adj. Factor Pilot            = ',rate_adjust
 write(*,'(A,1x,f4.2,A,f4.2)') ' Targeted Rej./Acc. rates              = ',acc_inf,'-',acc_sup
 write(*,'(A,1x,i6)')          ' R.N.G. seed                            = ',seed 

!!on rentre dans le dur

 if(opt_pheno) then
  allocate(PHENO_VAL(npop,npheno),dum_pheno_std(0:(nstream-1),npop))
  allocate(INITS_Delta(nmrk,npheno),INITS_Beta_i(nmrk,npheno),INITS_tau_beta(npheno),INIT_Pdelta(npheno))
  if(up_betai_coop) then
   allocate(acc_betai(nmrk,npheno),delta_betai(nmrk,npheno))
   delta_betai=0.05 ; acc_betai=0.
  end if
  if(opt_predcov) then
   allocate(PHENO_OBS(npop,npheno,2))
  end if
 end if

 if(opt_reg) then !forcement opt_pheno=TRUE donc pas besoin d'allouer
  allocate(PHENO_VAL(npop,npheno),dum_pheno_std(0:(nstream-1),npop))
  allocate(mean_rho_coef(nmrk,npheno,2),mean_isbf(nmrk,npheno,4),grid_is(0:(nstream-1),ninter_beta,npheno+1),dum_vect_grid(0:(nstream-1),ninter_beta))
  mean_rho_coef=0. !pearson
  mean_isbf=0. !1=m_bf ; 2=sd_bf ; 3=m_beta ; 4=sd_beta ! par importance sampling sur grille
  dum_real=(max_beta - min_beta)/(ninter_beta-1.)
  grid_is(:,1,1)=min_beta
  do dum_int=2,ninter_beta
   grid_is(:,dum_int,1)=grid_is(:,dum_int-1,1) + dum_real
  end do
 end if

 if(opt_pheno .or. opt_reg) then
  print *,''
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'        READING COVARIATE DATA'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,''
  if(opt_predcov) then ! on traite separement les deux cas (la deuxieme partie correspond exactement à l'ancienne version
   if(scale_cov) open(10,file=sum_phenostd_file,status='unknown') 
   open(1,file=Pheno_SDfile,status='old')
   open(2,file=Pheno_file,status='old') 
   do pheno=1,npheno
    read(2,*,iostat=err) PHENO_OBS(1:npop,pheno,1) 
    read(1,*,iostat=err) PHENO_OBS(1:npop,pheno,2)
    dum_int=npop
    do pop=1,npop
     if(PHENO_OBS(pop,pheno,2)<0) then
       PHENO_OBS(pop,pheno,1)=0. ; dum_int=dum_int-1
     end if
    end do
    dum_real = sum(PHENO_OBS(:,pheno,1))/real(dum_int)
    dum_real2 = sqrt((sum(PHENO_OBS(:,pheno,1)**2)/real(dum_int) - dum_real**2)*real(dum_int)/(real(dum_int)-1.))
    write(*,'(A,1x,i3,1x,A,f12.4,1x,A,f12.4,A)') 'Cov:',pheno,' (Mean =',dum_real,' ; SD =',dum_real2,')'
    if(scale_cov) then
     PHENO_OBS(:,pheno,1)=(PHENO_OBS(:,pheno,1)-dum_real )/dum_real2
     do pop=1,npop
      if(PHENO_OBS(pop,pheno,2)<0) PHENO_OBS(pop,pheno,1)=0. 
     end do
     write(10,'(1000(f10.6,1x))') PHENO_OBS(:,pheno,1)
    end if
    !transformation des sd en precision : permet de facilement tenir compte des pheno missing en mettant la precision a zero (le pheno etant deja a zero)
    do pop=1,npop
     if(PHENO_OBS(pop,pheno,2)<0) then
      PHENO_OBS(pop,pheno,2)=0.
     else
      PHENO_OBS(pop,pheno,2)=1./(PHENO_OBS(pop,pheno,2)*PHENO_OBS(pop,pheno,2))
     end if
    end do
   end do
    close(1) ; close(2) ; close(10)
  else
   open(1,file=Pheno_file,status='old') 
   if(scale_cov) open(10,file=sum_phenostd_file,status='unknown') 
   do pheno=1,npheno 
    read(1,*,iostat=err) PHENO_VAL(1:npop,pheno) 
    dum_real = sum(PHENO_VAL(:,pheno))/real(npop)
    dum_real2 = sqrt((sum(PHENO_VAL(:,pheno)**2)/npop - dum_real**2)*real(npop)/(real(npop)-1.))
    write(*,'(A,1x,i3,1x,A,f12.4,1x,A,f12.4,A)') 'Cov:',pheno,' (Mean =',dum_real,' ; SD =',dum_real2,')'
    if(scale_cov) then
     PHENO_VAL(:,pheno)=(PHENO_VAL(:,pheno)-dum_real )/dum_real2
     write(10,'(1000(f10.6,1x))') PHENO_VAL(:,pheno)
    end if
   end do
   close(1)
   if(scale_cov) close(10)
  end if
 end if

 if(opt_contrast) then
  print *,''
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'        READING CONTRASTS'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,''
  allocate(contrastes(ncontrast,npop),mean_contrast(ncontrast,nmrk,8),contrast_sd(ncontrast))
  mean_contrast=0.
  open(1,file=contrast_file,status='old') 
  do contrast=1,ncontrast 
   read(1,*,iostat=err) contrastes(contrast,1:npop)
   dum_int=0 ; dum_int2=0
   do pop=1,npop !verification
    if(contrastes(contrast,pop)/=-1 .and. contrastes(contrast,pop)/=0 .and. contrastes(contrast,pop)/=1) then
     print *,'ERROR: Only (integer) values of -1, 0 and 1 are allowed in the contrast file'
     stop
    else
     if(contrastes(contrast,pop)==1) dum_int=dum_int+1
     if(contrastes(contrast,pop)==-1) dum_int2=dum_int2+1
    end if
   end do
   if(dum_int==0 .or. dum_int2==0) then
    print *,'ERROR IN CONTRAST FILE for contrast: ',contrast
    print *,'At least one pop need to be 1 and one pop to be -1 for contrast to be computed'
    stop    
   else
    contrast_sd(contrast)=sqrt(real(dum_int)+real(dum_int2))
   end if
  end do
  close(1)
 end if
 
 !!!!!!!!!!!!!
 !!!Reading genotype file
 !!!!!!!!!!!!!
 print *,'' 
 if(indseq) then
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'        READING GENOTYPE LIKELIHOOD'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
  print *,''
  allocate(ind_pl(nmrk,npop,max_nind,3),cur_freq_like(nmrk,npop),dum_read_line(3*nind),dum_int_pop(npop),dum_real_pop(3)) 
  open(1,file=Geno_file,status='unknown') 
  read(1,'(a)',iostat=err) big_string
  do mrk=1,nmrk
   if(mod(mrk,10000)==0) print *,'Reading SNP',mrk
   dum_int_pop=0
   read(1,'(a)',iostat=err) big_string
  !BEGIN control of the last individual (control of all the individual too much time consuming)
   dum_int=0 ; dum_int2=0 !count blanks (must be = nind-1) and commas (must be = 2*nind)
   do ind=1,len_trim(big_string) 
    if(big_string(ind:ind)==' ') dum_int=dum_int+1
    if(big_string(ind:ind)==',') dum_int2=dum_int2+1
   end do
   if(dum_int+1/=nind) then
    write(*,'(A,i6)') 'ERROR parsing genotype likelihood file, line ',ind+1
    print *,trim(big_string)
    write(*,'(i6,A,i6,A)') dum_int,' spaces (separators of each individual genotype likelihood triplets) found instead of ',nind-1,' (nind-1) expected' 
    stop
   end if
   if(dum_int2 /= 2*nind) then
    write(*,'(A,i6)') 'ERROR parsing genotype likelihood file, line ',ind+1
    print *,trim(big_string)
    write(*,'(i6,A,i6,A)') dum_int2,' commas (separators of genotype likelihoods) found instead of ',2*nind,' (2*nind: 2 commas per individual genotype likelihood triplets) expected'
    write(*,'(A)') 'Reminder: only bi-allelic SNPs are supported by the model' 
    stop
   end if
   !END control
   read(big_string,*) dum_read_line !par dafaut comma et space=CS separator in fortran
   do ind=1,nind
    dum_real_pop=10.**(-0.1*dum_read_line((3*ind-2):(3*ind)))
    dum_real_pop=dum_real_pop/sum(dum_real_pop)
    dum_int_pop(ind_pop_index(ind))=dum_int_pop(ind_pop_index(ind)) + 1    
    ind_pl(mrk,ind_pop_index(ind),dum_int_pop(ind_pop_index(ind)),:)=dum_real_pop
   end do
  !  print *,i,1,pl(i,1,1,:),pl(i,39,8,:)!,trim(line)
  end do
  deallocate(dum_read_line,dum_int_pop,dum_real_pop)
  close(1)
 else
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'        READING COUNT DATA'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,''
  allocate(Y_OBS(nmrk,npop), N_OBS(nmrk,npop),missing_data(nmrk,npop),log_cn(nmrk,npop),dum_read_line(2*npop),dum_int_pop(npop))
  open(1,file=Geno_file,status='old') 
  do pop=1,npop
   dum_int_pop(pop)=2*pop-1 !on stocke la position
  end do
  do mrk=1,nmrk
   read(1,*,iostat=err) dum_read_line(1:(2*npop))
   Y_OBS(mrk,1:npop) = dum_read_line(dum_int_pop)
   N_OBS(mrk,1:npop) = Y_OBS(mrk,1:npop)  + dum_read_line(dum_int_pop+1) 
  end do
  deallocate(dum_read_line,dum_int_pop)
  close(1)
  write(*,'(A,1x,1000(i4,1x))') 'First SNP ref. Allele count ',Y_OBS(1,:) 
  write(*,'(A,1x,1000(i4,1x))') 'First SNP total Gene count  ',N_OBS(1,:) 
  write(*,'(A,1x,1000(i4,1x))') 'Last  SNP ref. Allele count ',Y_OBS(nmrk,:) 
  write(*,'(A,1x,1000(i4,1x))') 'Last  SNP total Gene count  ',N_OBS(nmrk,:)
!calcul des coef n=binomiaux utiles pour la suite (calcule de deviance...)
!a ce stade Y_obs peuvent etre des lectures ou des comptages
  missing_data=.false.
  do mrk=1,nmrk
   do pop=1,npop
    if(N_OBS(mrk,pop)<1) then
     missing_data(mrk,pop)=.true.
     log_cn(mrk,pop)=0.
    else 
     log_cn(mrk,pop)=log_binomial_coef(N_OBS(mrk,pop),Y_OBS(mrk,pop))
    end if
   end do
  end do
!!poolseq data
 if(poolseq) then !les Y_OBS et N_OBS deviennent des comptades et on introduit les lectures
  allocate(Y_READS(nmrk,npop),N_READS(nmrk,npop),delta_y(nmrk,npop),Yminmax(nmrk,npop,2),acc_y(nmrk,npop))
  Y_READS=Y_OBS ; N_READS=N_OBS ; delta_y=delta0_yij
  do pop=1,npop
   N_OBS(:,pop)=PoolSize(pop)
   do mrk=1,nmrk
     if(missing_data(mrk,pop)) then
       Yminmax(mrk,pop,1)=0 ; Yminmax(mrk,pop,2)=PoolSize(pop)
     else
      if(Y_READS(mrk,pop)/=0 .and. Y_READS(mrk,pop)/=N_READS(mrk,pop)) then
        Yminmax(mrk,pop,1)=1 ; Yminmax(mrk,pop,2)=PoolSize(pop)-1 !ymax=popsize-1 ; ymin=1 
      end if
      if(Y_READS(mrk,pop)==0) then
       Yminmax(mrk,pop,1)=0 ; Yminmax(mrk,pop,2)=PoolSize(pop)-1 !ymax=popsize-1 ; ymin=0
      end if
      if(Y_READS(mrk,pop)==N_READS(mrk,pop)) then
       Yminmax(mrk,pop,1)=1 ; Yminmax(mrk,pop,2)=PoolSize(pop)  ! ymax=popsize ; ymin=1
      end if
     end if
    end do
   end do
  end if
 end if

!!Allocation d'arrays necessaire pour la suite
 allocate(INITS_Pij(nmrk,npop,2) , INITS_Pi(nmrk) , INITS_LDA(npop,npop),omega_mat(npop,npop),&
         sum_covariates(nmrk,npop),delta_pi(nmrk),acc_pi(nmrk),muvar_cur(npop,2),&
         dum_real_pop(npop),dum_real_pop2(npop),c_mat(npop,npop),inv_c_mat(npop,npop),&
         array_real_pop(0:(nstream-1),npop,2),dum_freq(nmrk,npop,2))
 if(up_alpha_coop) then
   allocate(delta_pij(nmrk,1),acc_pij(nmrk,1)) 
 else
  allocate(delta_pij(nmrk,npop),acc_pij(nmrk,npop),omega_index_subpop(npop,npop-1),reg_mat(npop,npop-1),dum_int_pop(npop))
  dum_int_pop=(/(tst,tst=1,npop)/)
  do pop=1,npop
   omega_index_subpop(pop,:)=pack(dum_int_pop,dum_int_pop/=pop)
  end do
  deallocate(dum_int_pop)
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Info POP RACINE (le cas echeant) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(root_geno) then
  print *,''
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'        READING COUNT DATA for the root population'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,''
  allocate(root_missing(nmrk))
  open(1,file=Root_Geno_file,status='old')
  root_missing=.false.
  do mrk=1,nmrk
   read(1,*,iostat=err) R_CNT(mrk,1),R_CNT(mrk,2)
   R_CNT(mrk,3)=R_CNT(mrk,1)+R_CNT(mrk,2)
   if(R_CNT(mrk,3)==0) root_missing(mrk)=.true.
  end do
  close(1)
  if(root_pool) then
   allocate(R_READS(nmrk,3),delta_rooty(nmrk),RootYminmax(nmrk,2),acc_rooty(nmrk))
   R_READS=R_CNT ; delta_rooty=delta0_yij
   R_CNT(:,3)=root_poolsize
   do mrk=1,nmrk
    if(root_missing(mrk)) then
      RootYminmax(mrk,1)=0 ; RootYminmax(mrk,2)=root_poolsize
    else
     if(R_READS(mrk,1)/=0 .and. R_READS(mrk,2)/=0) then
       RootYminmax(mrk,1)=1 ; RootYminmax(mrk,2)=root_poolsize-1 !ymax=popsize-1 ; ymin=1 
     end if
     if(R_READS(mrk,1)==0) then
      RootYminmax(mrk,1)=0 ; RootYminmax(mrk,2)=root_poolsize-1 !ymax=popsize-1 ; ymin=0
     end if
     if(R_READS(mrk,2)==0) then
      RootYminmax(mrk,1)=1 ; RootYminmax(mrk,2)=root_poolsize  ! ymax=popsize ; ymin=1
     end if
    end if
   end do
  end if
 end if


 if(.not. estim_omega) then 
  print *,''
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'        READING OMEGA MATRIX'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,''
  open(1,file=omega_file,status='old') 
  do pop=1,npop
   read(1,*,iostat=err)  omega_mat(pop,1:npop) 
  end do
  close(1)
  INITS_LDA=inv(omega_mat)
  c_mat = omega_mat ; call chol(c_mat) 
  inv_c_mat=inv(c_mat)
  if(.not. up_alpha_coop) call compute_reg_mat()
 end if

 delta_pi=delta0_pi ; delta_pij=delta0_pij

end subroutine void_read_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!PRINT HELP PAGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_print_help()
   print *,"Version 2.2"
   print *,''
   print *,'Usage: BayPass [options]'
   print *,''
   print *,'Options:'
   print *,' I)   General Options:'
   write(*,'(A)') '  -help                   Display the help page '
   write(*,'(A)') '  -gfile           CHAR   Genotyping Data File                                            (always required)'
   write(*,'(A)') '  -efile           CHAR   Covariate file: activate Covariate Mode                         (def="")'
   write(*,'(A)') '  -scalecov        CHAR   Scale covariates                                                (def="")'
   write(*,'(A)') '  -contrastfile    CHAR   Contrast to be computed                                         (def="")'
   write(*,'(A)') '  -poolsizefile    CHAR   Name of the Pool Size file => activate PoolSeq mode             (def="")'
   write(*,'(A)') '  -outprefix       CHAR   Prefix used for the output files                                (def="")' 
   print *,''
   print *,' II)  Model Options:'
   write(*,'(A)') '  -omegafile       CHAR   Omega matrix file => inactivate estim. of omega                 (def="")'
   write(*,'(A)') '  -rho             INT    Rho parameter of the Wishart prior on omega                     (def=1)'
   write(*,'(A)') '  -setpibetapar           Inactivate estimation of the Pi beta priors parameters             '
   write(*,'(A)') '  -betapiprior     FLOAT2 Pi Beta prior parameters (if -setpibetapar)                     (def=1.0 1.0)'   
   write(*,'(A)') '  -minbeta         FLOAT  Lower beta coef. for the grid                                   (def=-0.3) '
   write(*,'(A)') '  -maxbeta         FLOAT  Upper beta coef. for the grid                                   (def= 0.3) '
   print *,''
   print *,'  I.1)  IS covariate mode (default covariate mode):'
   write(*,'(A)') '    -nbetagrid       INT    Number of grid points (IS covariate mode)                     (def=201) '
   print *,''
   print *,'  I.2)  MCMC covariate mode:'
   write(*,'(A)') '    -covmcmc                Activate mcmc covariate mode (desactivate estim. of omega) '
   write(*,'(A)') '    -auxmodel               Activate Auxiliary variable mode to estimate BF        '
   write(*,'(A)') '    -isingbeta       FLOAT  Beta (so-called inverse temperature) of the Ising model       (def=0.0)       '   
   write(*,'(A)') '    -auxPbetaprior   FLOAT2 auxiliary P Beta prior parameters                             (def=0.02 1.98)                         '   
   print *,''
   print *,' III) MCMC Options:'
   write(*,'(A)') '  -nthreads        INT    Number of threads                                               (def=1) '
   write(*,'(A)') '  -nval            INT    Number of post-burnin and thinned samples to generate           (def=1000) '
   write(*,'(A)') '  -thin            INT    Size of the thinning (record one every thin post-burnin sample) (def=20) '
   write(*,'(A)') '  -burnin          INT    Burn-in length                                                  (def=5000) '
   write(*,'(A)') '  -npilot          INT    Number of pilot runs (to adjust proposal distributions)         (def=20) '
   write(*,'(A)') '  -pilotlength     INT    Pilot run length                                                (def=500) '
   write(*,'(A)') '  -accinf          FLOAT  Lower target acceptance rate bound                              (def=0.25) '
   write(*,'(A)') '  -accsup          FLOAT  Upper target acceptance rate bound                              (def=0.40) '
   write(*,'(A)') '  -adjrate         FLOAT  Adjustement factor                                              (def=1.25) '
   write(*,'(A)') '  -d0pi            FLOAT  Initial delta for the pi all. freq. proposal                    (def=0.5)'
   write(*,'(A)') '  -upalphaalt             Alternative update of the pij  '
   write(*,'(A)') '  -d0pij           FLOAT  Initial delta for the pij all. freq. proposal (alt. update)     (def=0.05)'
   write(*,'(A)') '  -d0yij           INT    Initial delta for the yij all. count (PoolSeq mode)             (def=1)'
   write(*,'(A)') '  -d0cj            FLOAT  If nicholsonprior is set for Omega, initial delta for the cj    (def=0.05)'
   write(*,'(A)') '  -seed            INT    Random Number Generator seed                                    (def=5001) '
   write(*,'(A)') '  -print_omega_samples    Print post-burnin and thinned samples of Omega in a file ' 


end subroutine void_print_help

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!INITIALISATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_initialisation()
 !1) initialisation des Pij, Alpha_i et Pi_i
 INITS_Pij=0. ; INITS_Pi=0.0 !on initialise les INITS_Pi à la moyenne des A_IJ (sauf si info sur la racine)
 do mrk=1,nmrk
  do pop=1,npop
   if(indseq) then
    !A uniform prior is assumed on the three possible genotypes
    !Reference allele frequency
    do ind=1,nind_par_pop(pop)
     INITS_Pij(mrk,pop,1)=INITS_Pij(mrk,pop,1) + ind_pl(mrk,pop,ind,1) + 0.5*ind_pl(mrk,pop,ind,2)
    end do
    INITS_Pij(mrk,pop,1)=INITS_Pij(mrk,pop,1)/nind_par_pop(pop)
    !initialisation des likelihood
    !hwe_prob
    dum_real=max(1e-8,min(1.-1e-8,INITS_Pij(mrk,pop,1))) !Attention: eps doit matcher avec celui de update_alpha_indseq
    dum_real5(1)=dum_real*dum_real ; dum_real5(2)=2.*(1.-dum_real)*dum_real ; dum_real5(3)=(1.-dum_real)*(1.-dum_real)
    cur_freq_like(mrk,pop)=0.
    do ind=1,nind_par_pop(pop)
     cur_freq_like(mrk,pop)=cur_freq_like(mrk,pop)+log(dot_product(ind_pl(mrk,pop,ind,:),dum_real5(1:3)))
    end do    
   else
    if(poolseq) then !Cas du poolseq: initialisation des Y_OBS
     if(missing_data(mrk,pop)) then
      Y_OBS(mrk,pop)=PoolSize(pop)/2
     else
      dum_real=real(Y_READS(mrk,pop))*real(PoolSize(pop))/real(N_READS(mrk,pop))
      Y_OBS(mrk,pop)=nint(dum_real)
      if(Y_READS(mrk,pop)>0 .and. Y_OBS(mrk,pop)==0) Y_OBS(mrk,pop)=1
      if(Y_READS(mrk,pop)<N_READS(mrk,pop) .and. Y_OBS(mrk,pop)==PoolSize(pop)) Y_OBS(mrk,pop)=PoolSize(pop)-1
     end if
    end if
    INITS_Pij(mrk,pop,1)=(Y_OBS(mrk,pop)+1.)/(N_OBS(mrk,pop)+2.)
   end if
  end do
  INITS_Pi(mrk)=min(0.999,max(0.001,sum(INITS_Pij(mrk,:,1))/real(npop)))
  if(root_pool .and. (.not. root_missing(mrk))) then
  !dans le cas de donnees cnt sur la pop root n garde l'initisalisation de Pi comme avant (de toute facon pas de probleme aux brnes puisque toujours >0)
  !juste besoin de regarder pour les pools car il faut initialiser correctement les comptages
  !si pool et missing data: on update pas les r_cnt qui resteront toujours a zero
    R_CNT(mrk,1)=nint(INITS_Pi(mrk)*root_poolsize)
    if(R_READS(mrk,1)>0 .and. R_CNT(mrk,1)==0) R_CNT(mrk,1)=1
    if(R_READS(mrk,2)>0 .and. R_CNT(mrk,1)==root_poolsize) R_CNT(mrk,1)=root_poolsize-1     
    R_CNT(mrk,2)=root_poolsize-R_CNT(mrk,1) !Rq: R_CNT(,3) deja initialise à root_poolsize à la lecture
  end if
  INITS_Pij(mrk,:,2)= (INITS_Pij(mrk,:,1) - INITS_Pi(mrk))/sqrt(INITS_Pi(mrk)*(1.-INITS_Pi(mrk)))
 end do

 if(estim_omega) then
 !2) initialisation de LDA: comme un phylogenie en etoile cad diagonale=1/cj
  INITS_LDA=0. ; omega_mat=0. ; c_mat=0. ; inv_c_mat=0.
  do pop=1,npop
    do mrk=1,nmrk
     INITS_LDA(pop,pop)= INITS_LDA(pop,pop) + INITS_Pij(mrk,pop,2)**2 
    end do
    INITS_LDA(pop,pop) = real(nmrk)/INITS_LDA(pop,pop)
    omega_mat(pop,pop) = 1./INITS_LDA(pop,pop)
    c_mat(pop,pop)=sqrt(omega_mat(pop,pop))
    inv_c_mat(pop,pop)=1./c_mat(pop,pop)
  end do  
!~   c_mat = omega_mat ; call chol(c_mat) ; inv_c_mat=inv(c_mat)
  if(.not. up_alpha_coop) reg_mat=0. ! call compute_reg_mat()
 end if

!initialisation des parametres de la beta_pi
 if(up_params_beta) then
  tmp_a=sum(INITS_Pi)/real(nmrk) ; tmp_b=sum(INITS_Pi*INITS_Pi-tmp_a**2)/(real(nmrk)-1.)
  dum_real=(tmp_a*(1-tmp_a))/tmp_b - 1.
  pi_beta_params(1)=max(0.001,tmp_a*dum_real)  ; pi_beta_params(2)=max(0.001,(1.-tmp_a)*dum_real)  !estimateur des moments
  tmp_b=sum(pi_beta_params)
  pi_beta_update_save(3)=gamma_log(tmp_b)
  pi_beta_update_save(4)=gamma_log(tmp_b*tmp_a) ; pi_beta_update_save(5)=gamma_log(tmp_b*(1-tmp_a))
 end if
 
 sum_covariates=0.
 if(opt_pheno) then
  INITS_Delta = 1 ; INITS_Beta_i=0. ; INIT_Pdelta=P_beta_params(1)/sum(P_beta_params) ; INITS_tau_beta = tau_beta0
!~   if(opt_aux) INITS_tau_beta=INITS_tau_beta*INIT_Pdelta !=>on suppose que la transformation est deja dans la valeur passée en commande
  if(opt_predcov) PHENO_VAL=PHENO_OBS(:,:,1) 
 end if
 
end subroutine void_initialisation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!PRINT node values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_printval()
 if(poolseq) then
  tmp_min=minval(real(Y_OBS)) ; tmp_max=maxval(real(Y_OBS))
  write(*,'(A,(1x,f10.1,2x),A,(1x,f10.1,2x),A,(1x,f10.1))') 'Yij min         =',tmp_min,'Yij max         =',tmp_max,&
                                                            'Average Yij     =',sum(real(Y_OBS))/(real(nmrk)*real(npop))
 end if
 tmp_min=minval(INITS_Pij(:,:,1)) ; tmp_max=maxval(INITS_Pij(:,:,1))
 write(*,'(A,(1x,f12.5,2x),A,(1x,f12.5,2x),A,(1x,f12.5))') 'Pij min         =',tmp_min,'Pij max         =',tmp_max,&
                                                           'Average Pij     =',sum(INITS_Pij(:,:,1))/(real(nmrk)*real(npop))
 tmp_min=minval(INITS_Pij(:,:,2)) ; tmp_max=maxval(INITS_Pij(:,:,2))
 write(*,'(A,(1x,f12.5,2x),A,(1x,f12.5,2x),A,(1x,f12.5))') 'Pij_til min     =',tmp_min,'Pij_til max     =',tmp_max,&
                                                           'Average Pij_til =',sum(INITS_Pij(:,:,2))/(real(nmrk)*real(npop))
 tmp_min=minval(INITS_Pi) ; tmp_max=maxval(INITS_Pi)
 write(*,'(A,(1x,f12.5,2x),A,(1x,f12.5,2x),A,(1x,f12.5))') 'PI min          =',tmp_min,'PI max          =',tmp_max,&
                                                           'Average PI      =',sum(INITS_Pi)/real(nmrk)

 if(up_params_beta) write(*,'(A,1x,2(f12.5,1x))')          'Pi_beta_params values ',pi_beta_params

 if(printintmat) then
  print *,'OMEGA node values '
  do pop=1,npop
   write(*,'(1000(f12.8,1x))') ( omega_mat(pop,dum_int), dum_int=1,npop )
  end do
 end if

end subroutine void_printval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!AJUSTEMENT DES PROPOSALS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_adjust_proposals()
 tst_pij=0 ; tst_pi=0 ; tst_y=0 ; tst_b=0 ; tst_cj=0

 if(poolseq) then  !ajustement des y_ij
  do mrk=1,nmrk
   do pop=1,npop
    if(acc_y(mrk,pop)/pilot_length>acc_sup) then
!~      delta_y(mrk,pop)=nint(min(PoolSize(pop)+0.,delta_y(mrk,pop)+1.)) ; tst_y=tst_y+1
     delta_y(mrk,pop)=min(PoolSize(pop),delta_y(mrk,pop)+1) ; tst_y=tst_y+1
    end if
    if(acc_y(mrk,pop)/pilot_length<acc_inf) then
!~      delta_y(mrk,pop)=nint(max(1.,delta_y(mrk,pop)-1.)) ; tst_y=tst_y+1
     delta_y(mrk,pop)=max(1,delta_y(mrk,pop)-1) ; tst_y=tst_y+1
    end if
   end do
  end do
 write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_y,' out of ',nmrk*npop ,' y'
 end if

 if(root_pool) then
  dum_int=0 ; tst_y=0
  do mrk=1,nmrk
   if(.not. root_missing(mrk)) then
    dum_int=dum_int+1
    if(acc_rooty(mrk)/pilot_length>acc_sup) then
!~      delta_rooty(mrk)=nint(min(root_poolsize+0.,delta_rooty(mrk)+1.)) ; tst_y=tst_y+1
     delta_rooty(mrk)=min(root_poolsize,delta_rooty(mrk)+1) ; tst_y=tst_y+1
    end if
    if(acc_rooty(mrk)/pilot_length<acc_inf) then
!~      delta_rooty(mrk)=nint(max(1.,delta_rooty(mrk)-1.)) ; tst_y=tst_y+1
     delta_rooty(mrk)=max(1,delta_rooty(mrk)-1) ; tst_y=tst_y+1
    end if
   end if 
  end do
 write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_y,' out of ',dum_int ,' non null Root y'
 end if

 if(nich_prior) then
  do pop=1,npop
   if(acc_nichc(pop)/pilot_length>acc_sup) then
     delta_nichc(pop)=delta_nichc(pop)*rate_adjust ; tst_cj=tst_cj+1
   end if
   if(acc_nichc(pop)/pilot_length<acc_inf) then
     delta_nichc(pop)=delta_nichc(pop)/rate_adjust ; tst_cj=tst_cj+1
   end if
  end do 
  write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_cj,' out of ',npop ,' cj'  
 end if

!ajustement des alpha proposals
 dum_int=npop
 if(up_alpha_coop) dum_int=1
 do mrk=1,nmrk
  do pop=1,dum_int
   if(acc_pij(mrk,pop)/pilot_length>acc_sup) then
     delta_pij(mrk,pop)=delta_pij(mrk,pop)*rate_adjust ; tst_pij=tst_pij+1
   end if
   if(acc_pij(mrk,pop)/pilot_length<acc_inf) then
    delta_pij(mrk,pop)=delta_pij(mrk,pop)/rate_adjust ; tst_pij=tst_pij+1
   end if
  end do
 end do
 if(up_alpha_coop) then
  write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_pij,' out of ',nmrk ,' p_ij vectors'
 else
  write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_pij,' out of ',nmrk*npop ,' p_ij'
 end if

!ajustement des pi proposals (le cas echeant)
! if(.not. up_slc_pi) then
  do mrk=1,nmrk
   if(acc_pi(mrk)/pilot_length>acc_sup) then
    delta_pi(mrk)=delta_pi(mrk)*rate_adjust ; tst_pi=tst_pi+1
   end if    
   if(acc_pi(mrk)/pilot_length<acc_inf) then
    delta_pi(mrk)=delta_pi(mrk)/rate_adjust ; tst_pi=tst_pi+1
   end if  
  end do
  write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_pi,' out of ',nmrk ,' pi'
! end if

 if(opt_pheno .and. up_betai_coop) then
  do mrk=1,nmrk
   do pheno=1,npheno
    if(acc_betai(mrk,pheno)/pilot_length>acc_sup) then
     delta_betai(mrk,pheno)=delta_betai(mrk,pheno)*rate_adjust ; tst_b=tst_b+1
    end if    
    if(acc_betai(mrk,pheno)/pilot_length<acc_inf) then
     delta_betai(mrk,pheno)=delta_betai(mrk,pheno)/rate_adjust ; tst_b=tst_b+1
    end if  
   end do
  end do
  write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_b,' out of ',nmrk*npheno ,' betai'
 end if


 if(poolseq) write(*,'(A,(1x,f7.5,2x),A,(1x,f6.1))') '      Mean Acceptance Rate Y  = ',sum(acc_y)/(real(pilot_length)*real(nmrk)*real(npop)),&
                                             ' mean delta_y= ',sum(delta_y+0.)/(real(nmrk)*real(npop))
 if(up_alpha_coop) then
  write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Pij  = ',sum(acc_pij)/(real(pilot_length)*real(nmrk)),&
                                         ' mean delta_p= ',sum(delta_pij)/real(nmrk)
 else
  write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Pij  = ',sum(acc_pij)/(real(pilot_length)*real(nmrk)*real(npop)),&
                                         ' mean delta_p= ',sum(delta_pij)/(real(nmrk)*real(npop))
 end if
! if(.not. up_slc_pi) write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Pi   = ',sum(acc_pi)/(pilot_length*nmrk),&
!                                         ' mean delta_pi= ',sum(delta_pi)/nmrk
 write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Pi   = ',sum(acc_pi)/(real(pilot_length)*real(nmrk)),&
                                         ' mean delta_pi= ',sum(delta_pi)/nmrk

 if(opt_pheno .and. up_betai_coop) then
   write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Betai   = ',sum(acc_betai)/(real(pilot_length)*real(nmrk)*real(npheno)),&
                                          ' mean delta_betai= ',sum(delta_betai)/(real(nmrk)*real(npheno))
 end if


 if(up_params_beta) then !ajustement des beta_pi_mu et beta_pi_phi
  if(acc_beta_params(1)/pilot_length>acc_sup) then
     delta_pi_beta_mu=delta_pi_beta_mu*rate_adjust
  end if    
  if(acc_beta_params(1)/pilot_length<acc_inf) then
     delta_pi_beta_mu=delta_pi_beta_mu/rate_adjust
  end if 
 write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate b_mu = ',acc_beta_params(1)/real(pilot_length),&
                                         '      delta  = ',delta_pi_beta_mu
  if(acc_beta_params(2)/pilot_length>acc_sup) then
    delta_pi_beta_phi=delta_pi_beta_phi*rate_adjust
  end if    
  if(acc_beta_params(2)/pilot_length<acc_inf) then
    delta_pi_beta_phi=delta_pi_beta_phi/rate_adjust
  end if   
 write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate b_phi= ',acc_beta_params(2)/real(pilot_length),&
                                        '      delta  = ',delta_pi_beta_phi
 end if

end subroutine void_adjust_proposals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!MCMC_ITER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_mcmc_iter()
 if(poolseq) then
 !$OMP PARALLEL private(pop,cur_stream,dum_int,y_out)
 !$OMP do schedule(guided) 
  do mrk=1,nmrk
   cur_stream=omp_get_thread_num()
   do pop=1,npop
    if(missing_data(mrk,pop)) then
     acc_y(mrk,pop)=acc_y(mrk,pop)+1.
     if(INITS_Pij(mrk,pop,1)>0. .and. INITS_Pij(mrk,pop,1)<1.) then
      dum_int=random_binomial2(PoolSize(pop),INITS_Pij(mrk,pop,1),mts(cur_stream))
      Y_OBS(mrk,pop)=dum_int
      if(Y_OBS(mrk,pop)<0) print *,mrk,pop,Y_OBS(mrk,pop),missing_data(mrk,pop),INITS_Pij(mrk,pop,1),PoolSize(pop),dum_int 
     else
      if(INITS_Pij(mrk,pop,1)<1e-16) Y_OBS(mrk,pop)=0
      if(INITS_Pij(mrk,pop,1)>0.9999) Y_OBS(mrk,pop)=PoolSize(pop)    
     end if
    else
     call update_y(mts(cur_stream),Y_OBS(mrk,pop),PoolSize(pop),Y_READS(mrk,pop),N_READS(mrk,pop),&
                   Yminmax(mrk,pop,:),max(0.,min(1.,INITS_Pij(mrk,pop,1))),delta_y(mrk,pop),dum_int,y_out)
     acc_y(mrk,pop)=acc_y(mrk,pop)+dum_int ; Y_OBS(mrk,pop)=y_out
    end if
!~     if(Y_OBS(mrk,pop)<0) print *,mrk,pop,Y_OBS(mrk,pop),missing_data(mrk,pop),INITS_Pij(mrk,pop,1),PoolSize(pop),random_binomial2(PoolSize(pop),INITS_Pij(mrk,pop,1),mts(cur_stream)) 
   end do
  end do
  !$OMP END do   
  !$OMP END PARALLEL  
 end if

 if(root_pool) then
 !$OMP PARALLEL private(cur_stream,dum_int,y_out)
 !$OMP do schedule(guided) 
  do mrk=1,nmrk
   cur_stream=omp_get_thread_num()
    if(.not. root_missing(mrk)) then !si missing on update pas et R_CNT reste à zero
     call update_y(mts(cur_stream),R_CNT(mrk,1),root_poolsize,R_READS(mrk,1),R_READS(mrk,3),&
                   RootYminmax(mrk,:),INITS_Pi(mrk),delta_rooty(mrk),dum_int,y_out)
     acc_rooty(mrk)=acc_rooty(mrk)+dum_int ; R_CNT(mrk,1)=y_out ; R_CNT(mrk,2)=root_poolsize-y_out
    end if
  end do
  !$OMP END do   
  !$OMP END PARALLEL  
 end if

  if(up_alpha_coop) then
 !$OMP PARALLEL private(cur_stream,dum_int)
 !$OMP do schedule(guided) 
   do mrk=1,nmrk
    cur_stream=omp_get_thread_num()
!    print *,cur_stream,size(dum_real_pop2) => illustre probleme avec ifort et -openmp: siz(dum_real_pop2)=0!!!
    call update_alpha_vect(mts(cur_stream),Y_OBS(mrk,:),N_OBS(mrk,:),INITS_Pij(mrk,:,1),INITS_Pij(mrk,:,2),INITS_Pi(mrk),&
                           sum_covariates(mrk,:),delta_pij(mrk,1),INITS_LDA,c_mat,dum_int,array_real_pop(cur_stream,:,1),array_real_pop(cur_stream,:,2))
     INITS_Pij(mrk,:,1)=array_real_pop(cur_stream,:,1) ; INITS_Pij(mrk,:,2)=array_real_pop(cur_stream,:,2) ; acc_pij(mrk,1)=acc_pij(mrk,1)+dum_int 
    end do
  !$OMP END do   
  !$OMP END PARALLEL  
  else
   if(indseq) then
    do pop=1,npop
     alpha_up_muvar_param(2)=omega_mat(pop,pop)-dot_product(reg_mat(pop,:),omega_mat(pop,omega_index_subpop(pop,:)))
 !$OMP PARALLEL private(cur_stream,dum_int,dum_real,dum_real2)
 !$OMP do schedule(guided) 
     do mrk=1,nmrk
      cur_stream=omp_get_thread_num()
      alpha_up_muvar_param(1)=dot_product(reg_mat(pop,:),INITS_Pij(mrk,omega_index_subpop(pop,:),2))
      call update_alpha_indseq (mts(cur_stream),ind_pl(mrk,pop,1:nind_par_pop(pop),:),cur_freq_like(mrk,pop),INITS_Pij(mrk,pop,1),INITS_Pij(mrk,pop,2),INITS_Pi(mrk),&
                               sum_covariates(mrk,pop),alpha_up_muvar_param(1),alpha_up_muvar_param(2),delta_pij(mrk,pop),dum_int,dum_real,dum_real2,tmp_a)
      INITS_Pij(mrk,pop,1)=dum_real ; INITS_Pij(mrk,pop,2)=dum_real2 ; acc_pij(mrk,pop)=acc_pij(mrk,pop)+dum_int ; cur_freq_like(mrk,pop)=tmp_a
     end do
  !$OMP END do   
  !$OMP END PARALLEL 
    end do 
   else
    do pop=1,npop
     alpha_up_muvar_param(2)=omega_mat(pop,pop)-dot_product(reg_mat(pop,:),omega_mat(pop,omega_index_subpop(pop,:)))
 !$OMP PARALLEL private(cur_stream,dum_int,dum_real,dum_real2)
 !$OMP do schedule(guided) 
    do mrk=1,nmrk
     cur_stream=omp_get_thread_num()
     alpha_up_muvar_param(1)=dot_product(reg_mat(pop,:),INITS_Pij(mrk,omega_index_subpop(pop,:),2))
     call update_alpha (mts(cur_stream),Y_OBS(mrk,pop),N_OBS(mrk,pop),INITS_Pij(mrk,pop,1),INITS_Pij(mrk,pop,2),INITS_Pi(mrk),&
                        sum_covariates(mrk,pop),alpha_up_muvar_param(1),alpha_up_muvar_param(2),delta_pij(mrk,pop),dum_int,dum_real,dum_real2)
      INITS_Pij(mrk,pop,1)=dum_real ; INITS_Pij(mrk,pop,2)=dum_real2 ; acc_pij(mrk,pop)=acc_pij(mrk,pop)+dum_int 
     end do
  !$OMP END do   
  !$OMP END PARALLEL 
    end do
   end if
  end if

!  if(up_slc_pi) then
!   do mrk=1,nmrk
!    INITS_Pi(mrk)=slc_pi(INITS_Pi(mrk),INITS_Pij(mrk,:,:),INITS_LDA,pi_beta_params) 
!    INITS_Pij(mrk,:,2)= (INITS_Pij(mrk,:,1) - INITS_Pi(mrk))/sqrt(INITS_Pi(mrk)*(1-INITS_Pi(mrk)))
!   end do 
!  else
 !$OMP PARALLEL private(cur_stream,dum_int,dum_real)
 !$OMP do schedule(guided) 
   do mrk=1,nmrk
    cur_stream=omp_get_thread_num()
    call update_pi(mts(cur_stream),INITS_Pi(mrk),INITS_Pij(mrk,:,:),sum_covariates(mrk,:),INITS_LDA,R_CNT(mrk,1:2),pi_beta_params,delta_pi(mrk),dum_int,dum_real,array_real_pop(cur_stream,:,1))
    acc_pi(mrk)=acc_pi(mrk)+dum_int ; INITS_Pi(mrk)=dum_real ; INITS_Pij(mrk,:,2)= array_real_pop(cur_stream,:,1)
   end do 
  !$OMP END do   
  !$OMP END PARALLEL 

  if(up_params_beta) then
   cur_stream=0
   pi_beta_update_save(1:2)=0.
   do mrk=1,nmrk
     pi_beta_update_save(1)=pi_beta_update_save(1) + log(INITS_Pi(mrk))
     pi_beta_update_save(2)=pi_beta_update_save(2) + log(1.-INITS_Pi(mrk))
   end do
    call update_beta_pi_params(mts(cur_stream),pi_beta_params,nmrk,pi_beta_update_save,delta_pi_beta_mu,&
                                delta_pi_beta_phi,dum_int,dum_int2,dum_real5)
    acc_beta_params(1)=acc_beta_params(1) + dum_int ; acc_beta_params(2)=acc_beta_params(2) + dum_int2
    pi_beta_params=dum_real5(1:2) ; pi_beta_update_save(3:5)=dum_real5(3:5)
  end if

  if(estim_omega) then
   if(nich_prior) then
   !$OMP PARALLEL private(cur_stream,dum_int,dum_real)
   !$OMP do schedule(guided) 
    do pop=1,npop
     cur_stream=omp_get_thread_num()
     call update_nichc(mts(cur_stream),nmrk,omega_mat(pop,pop),sum(INITS_Pij(:,pop,2)**2),delta_nichc(pop),dum_int,dum_real)
     acc_nichc(pop)=acc_nichc(pop)+dum_int ; omega_mat(pop,pop)=dum_real ; INITS_LDA(pop,pop)=1./dum_real 
     c_mat(pop,pop)=sqrt(omega_mat(pop,pop)) ; inv_c_mat(pop,pop)=1./c_mat(pop,pop) !on touche pas a reg_mat qui reste a 0 (valeur initiale)
    end do 
   !$OMP END do   
   !$OMP END PARALLEL 
   else
    cur_stream=0
    INITS_LDA=up_lda_mat() ; omega_mat=inv(INITS_LDA)
    c_mat=omega_mat ; call chol(c_mat)
    inv_c_mat=inv(c_mat)
    if(.not. up_alpha_coop) call compute_reg_mat()
   end if
  end if

  if(opt_pheno) then
   do pheno=1,npheno
    if(up_betai_coop) then
 !$OMP PARALLEL private(cur_stream,dum_int,dum_real)
 !$OMP do schedule(guided) 
     do mrk=1,nmrk
      cur_stream=omp_get_thread_num()
      call up_beta_i_coop(mts(cur_stream),INITS_Beta_i(mrk,pheno),npop,INITS_Pij(mrk,:,1),inv_c_mat,INITS_Delta(mrk,pheno),PHENO_VAL(:,pheno),&
                          sum_covariates(mrk,:),INITS_Pi(mrk),min_beta,max_beta,delta_betai(mrk,pheno),dum_real,array_real_pop(cur_stream,:,1),dum_int) 
      INITS_Beta_i(mrk,pheno)=dum_real ; sum_covariates(mrk,:)=array_real_pop(cur_stream,:,1) ; acc_betai(mrk,pheno)=acc_betai(mrk,pheno) + dum_int
     end do
  !$OMP END do   
  !$OMP END PARALLEL 
    else
 !$OMP PARALLEL private(cur_stream,dum_real)
 !$OMP do schedule(guided) 
     do mrk=1,nmrk
      cur_stream=omp_get_thread_num()
      call up_beta_i(mts(cur_stream),INITS_Beta_i(mrk,pheno),npop,INITS_Pij(mrk,:,1),inv_c_mat,INITS_Delta(mrk,pheno),PHENO_VAL(:,pheno),&
                          sum_covariates(mrk,:),INITS_Pi(mrk),INITS_tau_beta(pheno),dum_real,array_real_pop(cur_stream,:,1)) 
      INITS_Beta_i(mrk,pheno)=dum_real ; sum_covariates(mrk,:)=array_real_pop(cur_stream,:,1)
     end do
  !$OMP END do   
  !$OMP END PARALLEL 
    end if
    if(opt_aux .and. out_pilot) then
     if(opt_ising) then
      cur_stream=0
      do mrk=1,nmrk
       call up_delta_ising(mts(cur_stream),mrk,INITS_Delta(:,pheno),ising_b,INITS_Beta_i(mrk,pheno),npop,nmrk,INIT_Pdelta(pheno),&
                          INITS_Pij(mrk,:,1),inv_c_mat,PHENO_VAL(:,pheno),sum_covariates(mrk,:),INITS_Pi(mrk),dum_int,dum_real_pop)  
       INITS_Delta(mrk,pheno)=dum_int ; sum_covariates(mrk,:)=dum_real_pop
      end do
     else
 !$OMP PARALLEL private(cur_stream,dum_int)
 !$OMP do schedule(guided) 
      do mrk=1,nmrk
       cur_stream=omp_get_thread_num()
       call up_delta(mts(cur_stream),INITS_Delta(mrk,pheno),INITS_Beta_i(mrk,pheno),npop,INIT_Pdelta(pheno),INITS_Pij(mrk,:,1),inv_c_mat,&
                     PHENO_VAL(:,pheno),sum_covariates(mrk,:),INITS_Pi(mrk),dum_int,array_real_pop(cur_stream,:,1))  
       INITS_Delta(mrk,pheno)=dum_int ; sum_covariates(mrk,:)=array_real_pop(cur_stream,:,1)
      end do
  !$OMP END do   
  !$OMP END PARALLEL 
     end if
    end if
   end do

    if(estim_tau_beta) then
     cur_stream=0
     do pheno=1,npheno
      INITS_tau_beta(pheno)=up_tau_beta(mts(cur_stream),INITS_Beta_i(:,pheno),INITS_Delta(:,pheno),k_tau_prior,l_tau_prior)
     end do
    end if
    
    if(opt_aux .and. out_pilot) then
     cur_stream=0 !pas la peine de paralleliser ca
     do pheno=1,npheno
      INIT_Pdelta(pheno)=up_P(mts(cur_stream),INITS_Delta(:,pheno),P_beta_params(1),P_beta_params(2))
     end do
    end if

    if(opt_predcov) then
      cur_stream=0
      do pheno=1,npheno !ne pas paralleliser car parallelisation interne a up_phenoz
       call up_phenoz(mts(cur_stream),PHENO_OBS(:,pheno,:),covgaussprior,nmrk,npop,PHENO_VAL(:,pheno),INITS_Pij(:,:,1),INITS_Delta(:,pheno)*INITS_Beta_i(:,pheno),&
                      INITS_Pi,inv_c_mat,sum_covariates,array_real_pop(cur_stream,:,1))
       PHENO_VAL(:,pheno)=array_real_pop(cur_stream,:,1)               
      end do
    end if
  end if

end subroutine void_mcmc_iter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! CALCUL SUMMARY STATS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_update_summary()

if(.not. indseq) then ! calcul de la devaince pour donnees individuelles a voir plus tard
if(poolseq) then
 do pop=1,npop
  dum_freq(:,pop,1)=Y_OBS(:,pop)*1./PoolSize(pop)
 end do
 mean_deviance=(mean_deviance*(real(iter)-1.) + deviance(Y_READS,N_READS,dum_freq(:,:,1),log_cn))/real(iter)
 !$OMP PARALLEL private(pop)
 !$OMP do schedule(guided) 
 do mrk=1,nmrk
  do pop=1,npop
   CPO(mrk,pop)=(CPO(mrk,pop)*(real(iter)-1.) + &
                exp(-1.*logl_ij_beta(Y_READS(mrk,pop),N_READS(mrk,pop),dum_freq(mrk,pop,1),log_cn(mrk,pop))))/real(iter)
  end do
 end do
!$OMP END do   
!$OMP END PARALLEL 
else
 mean_deviance=(mean_deviance*(real(iter)-1.) + deviance(Y_OBS,N_OBS,INITS_Pij(:,:,1),log_cn))/real(iter)
 !$OMP PARALLEL private(pop)
 !$OMP do schedule(guided) 
 do mrk=1,nmrk
  do pop=1,npop
   CPO(mrk,pop)=(CPO(mrk,pop)*(real(iter)-1.) + &
                exp(-1.*logl_ij_beta(Y_OBS(mrk,pop),N_OBS(mrk,pop),INITS_Pij(mrk,pop,1),log_cn(mrk,pop))))/real(iter)
  end do
 end do
!$OMP END do   
!$OMP END PARALLEL 
end if
end if

!$OMP PARALLEL private(pop,contrast,dum_real)
!$OMP do schedule(guided) 
do mrk=1,nmrk
! cur_stream=omp_get_thread_num()
 mean_pi(mrk,1)=(mean_pi(mrk,1)*(real(iter)-1.)+INITS_Pi(mrk))/real(iter)
 mean_pi(mrk,2)=(mean_pi(mrk,2)*(real(iter)-1.)+(INITS_Pi(mrk))**2)/real(iter)
 !pij standardisée et xtx stat
 dum_freq(mrk,:,1) = matmul(inv_c_mat,INITS_Pij(mrk,:,2))
 dum_real=dot_product(dum_freq(mrk,:,1),dum_freq(mrk,:,1))
! dum_real=dot_product(INITS_Pij(mrk,:,2),matmul(INITS_LDA,INITS_Pij(mrk,:,2)))
 mean_xtx(mrk,1)=(mean_xtx(mrk,1)*(real(iter)-1.)+dum_real)/real(iter)
 mean_xtx(mrk,2)=(mean_xtx(mrk,2)*(real(iter)-1.)+dum_real**2)/real(iter)
 do pop=1,npop
  if(poolseq) then
    mean_yij(mrk,pop,1)=(mean_yij(mrk,pop,1)*(real(iter)-1.)+real(Y_OBS(mrk,pop)))/real(iter)
    mean_yij(mrk,pop,2)=(mean_yij(mrk,pop,2)*(real(iter)-1.)+real(Y_OBS(mrk,pop))**2)/real(iter)
  end if
  dum_real=max(0.,min(1.,INITS_Pij(mrk,pop,1)))
!~   dum_freq(mrk,pop,2)=(dum_real-INITS_Pi(mrk))/sqrt(INITS_Pi(mrk)*(1.-INITS_Pi(mrk))) ! pour les contrastes 2 (sinon pas utile)
  mean_pij(mrk,pop,1)=(mean_pij(mrk,pop,1)*(real(iter)-1.)+dum_real)/real(iter)
  mean_pij(mrk,pop,2)=(mean_pij(mrk,pop,2)*(real(iter)-1.)+dum_real**2)/real(iter)
  mean_pij(mrk,pop,3)=(mean_pij(mrk,pop,3)*(real(iter)-1.)+dum_freq(mrk,pop,1))/real(iter)
  mean_pij(mrk,pop,4)=(mean_pij(mrk,pop,4)*(real(iter)-1.)+dum_freq(mrk,pop,1)**2)/real(iter)
 end do
 dum_freq(mrk,:,2) = matmul(inv_c_mat,dum_freq(mrk,:,2)) ! pour les contrastes 2 (sinon pas utile)
 if(root_pool) then
    mean_rooty(mrk,1)=(mean_rooty(mrk,1)*(real(iter)-1.)+real(R_CNT(mrk,1)))/real(iter)
    mean_rooty(mrk,2)=(mean_rooty(mrk,2)*(real(iter)-1.)+real(R_CNT(mrk,1))**2.)/real(iter) 
 end if
 if(opt_contrast) then
  do contrast=1,ncontrast
   dum_real=sum(dum_freq(mrk,:,1)*contrastes(contrast,:)) / contrast_sd(contrast)
   mean_contrast(contrast,mrk,1)=(mean_contrast(contrast,mrk,1)*(real(iter)-1.)+dum_real)/real(iter)
   mean_contrast(contrast,mrk,2)=(mean_contrast(contrast,mrk,2)*(real(iter)-1.)+dum_real**2)/real(iter)
   dum_real=dum_real**2
   mean_contrast(contrast,mrk,3)=(mean_contrast(contrast,mrk,3)*(real(iter)-1.)+dum_real)/real(iter)
   mean_contrast(contrast,mrk,4)=(mean_contrast(contrast,mrk,4)*(real(iter)-1.)+dum_real**2)/real(iter)   
!~    dum_real=sum(dum_freq(mrk,:,2)*contrastes(contrast,:)) / contrast_sd(contrast)
!~    mean_contrast(contrast,mrk,5)=(mean_contrast(contrast,mrk,5)*(real(iter)-1.)+dum_real)/real(iter)
!~    mean_contrast(contrast,mrk,6)=(mean_contrast(contrast,mrk,6)*(real(iter)-1.)+dum_real**2)/real(iter)
!~    dum_real=dum_real**2
!~    mean_contrast(contrast,mrk,7)=(mean_contrast(contrast,mrk,7)*(real(iter)-1.)+dum_real)/real(iter)
!~    mean_contrast(contrast,mrk,8)=(mean_contrast(contrast,mrk,8)*(real(iter)-1.)+dum_real**2)/real(iter)   
  end do
 end if 
end do
!$OMP END do   
!$OMP END PARALLEL 

!!!!!!!!!!!!!!!!!!!!!!!
!!calculs à partir des pij standardisées et SCALEES PAR POP: XtX*, etc.: ATTENTION: verifier que dum_freq n'a pas change!!!!
!sert à rien car tres peu de differences par rapport aux contrastes classique ou au XtX
!En effet le fait de moyenner recree le shrinkage
!~  do pop=1,npop
!~   dum_real=sum(dum_freq(:,pop))/nmrk
!~   dum_real2=sqrt(sum((dum_freq(:,pop)-dum_real)**2)/(nmrk-1.))
!~   dum_freq(:,pop)=(dum_freq(:,pop) - dum_real)/dum_real2
!~  end do
!~ dum_real=sum(dum_freq(:,:,1))/(real(nmrk)*real(npop))
!~ dum_real2=sqrt(sum((dum_freq(:,:,1)-dum_real)**2)/(real(nmrk)*real(npop)-1.))
!~ dum_freq(:,:,1)=(dum_freq(:,:,1) - dum_real)/dum_real2
!~ !$OMP PARALLEL private(pop,contrast,dum_real)
!~ !$OMP do schedule(guided) 
!~ do mrk=1,nmrk
!~  dum_real=dot_product(dum_freq(mrk,:,1),dum_freq(mrk,:,1))
!~  mean_xtx(mrk,3)=(mean_xtx(mrk,3)*(iter-1)+dum_real)/iter
!~  mean_xtx(mrk,4)=(mean_xtx(mrk,4)*(iter-1)+dum_real**2)/iter
!~  !!contrastes le cas echeant à partir des valeur scalees et non scalees
!~  if(opt_contrast) then
!~   do contrast=1,ncontrast
!~    dum_real=((sum(dum_freq(mrk,:,1)*contrastes(contrast,:)))) / contrast_sd(contrast)
!~    mean_contrast(contrast,mrk,5)=(mean_contrast(contrast,mrk,5)*(iter-1)+dum_real)/iter
!~    mean_contrast(contrast,mrk,6)=(mean_contrast(contrast,mrk,6)*(iter-1)+dum_real**2)/iter
!~   end do
!~  end if
!~ end do
!~ !$OMP END do   
!~ !$OMP END PARALLEL 


if(up_params_beta) then
 mean_beta_params(1,1)=(mean_beta_params(1,1)*(real(iter)-1.)+pi_beta_params(1))/real(iter)
 mean_beta_params(1,2)=(mean_beta_params(1,2)*(real(iter)-1.)+(pi_beta_params(1))**2)/real(iter)
 mean_beta_params(2,1)=(mean_beta_params(2,1)*(real(iter)-1.)+pi_beta_params(2))/real(iter)
 mean_beta_params(2,2)=(mean_beta_params(2,2)*(real(iter)-1.)+(pi_beta_params(2))**2)/real(iter)
end if

 mean_lda(:,:,1)=(mean_lda(:,:,1)*(real(iter)-1.)+ INITS_LDA)/real(iter)
 mean_lda(:,:,2)=(mean_lda(:,:,2)*(real(iter)-1.)+ INITS_LDA**2)/real(iter)
 mean_omega(:,:,1)=(mean_omega(:,:,1)*(real(iter)-1.)+ omega_mat)/real(iter)
 mean_omega(:,:,2)=(mean_omega(:,:,2)*(real(iter)-1.)+ omega_mat**2)/real(iter)

 if(opt_pheno) then
  if(opt_aux) then
   do pheno=1,npheno
    do mrk=1,nmrk
     dum_real=INITS_Beta_i(mrk,pheno)*INITS_Delta(mrk,pheno)
     mean_beta_i(mrk,pheno,1)=(mean_beta_i(mrk,pheno,1)*real(iter-1)+dum_real)/real(iter)
     mean_beta_i(mrk,pheno,2)=(mean_beta_i(mrk,pheno,2)*real(iter-1)+dum_real**2)/real(iter)
     mean_delta_i(mrk,pheno)=(mean_delta_i(mrk,pheno)*real(iter-1)+INITS_Delta(mrk,pheno))/real(iter)
    end do
   mean_p_delta(pheno,1)=(mean_p_delta(pheno,1)*real(iter-1)+INIT_Pdelta(pheno))/real(iter)
   mean_p_delta(pheno,2)=(mean_p_delta(pheno,2)*real(iter-1)+INIT_Pdelta(pheno)**2)/real(iter)
   end do
  else
   do pheno=1,npheno
    do mrk=1,nmrk
     mean_beta_i(mrk,pheno,1)=(mean_beta_i(mrk,pheno,1)*real(iter-1)+INITS_Beta_i(mrk,pheno))/real(iter)
     mean_beta_i(mrk,pheno,2)=(mean_beta_i(mrk,pheno,2)*real(iter-1)+INITS_Beta_i(mrk,pheno)**2)/real(iter)
    end do
    if(estim_tau_beta) then
     mean_tau_beta(pheno,1)=(mean_tau_beta(pheno,1)*real(iter-1)+ INITS_tau_beta(pheno))/real(iter)
     mean_tau_beta(pheno,2)=(mean_tau_beta(pheno,2)*real(iter-1)+ INITS_tau_beta(pheno)**2)/real(iter)
    end if 
   end do
  end if
  if(opt_predcov) then
   do pheno=1,npheno
     mean_phenoval(:,pheno,1)=(mean_phenoval(:,pheno,1)*real(iter-1)+ PHENO_VAL(:,pheno))/real(iter)
     mean_phenoval(:,pheno,2)=(mean_phenoval(:,pheno,2)*real(iter-1)+ PHENO_VAL(:,pheno)**2)/real(iter)
   end do
  end if
 end if

 if(opt_reg) then
  do pheno=1,npheno
   dum_real_pop2 = matmul(inv_c_mat,PHENO_VAL(:,pheno))
!$OMP PARALLEL private(cur_stream,dum_int,dum_real)
!$OMP do schedule(guided) 
   do mrk=1,nmrk
    cur_stream=omp_get_thread_num()
    dum_pheno_std(cur_stream,:)=dum_real_pop2 / sqrt(INITS_Pi(mrk)*(1-INITS_Pi(mrk)))
    do dum_int=1,ninter_beta
     grid_is(cur_stream,dum_int,pheno+1)= -0.5*(grid_is(cur_stream,dum_int,1)**2)*(sum(dum_pheno_std(cur_stream,:)**2))
    end do
    array_real_pop(cur_stream,:,1) = matmul(inv_c_mat,INITS_Pij(mrk,:,2))
    dum_real= rho_pearson(dum_pheno_std(cur_stream,:),array_real_pop(cur_stream,:,1)) ! si on veut faire comme dans G&C (2013): rho_pearson(dum_real_pop2,array_real_pop(cur_stream,:,1))
    mean_rho_coef(mrk,pheno,1)=(mean_rho_coef(mrk,pheno,1)*real(iter-1)+ dum_real)/real(iter)
    mean_rho_coef(mrk,pheno,2)=(mean_rho_coef(mrk,pheno,2)*real(iter-1)+ dum_real**2)/real(iter)
    !calcul des stats par importance
    do dum_int=1,ninter_beta
     dum_vect_grid(cur_stream,dum_int)=exp(grid_is(cur_stream,dum_int,1)*sum(array_real_pop(cur_stream,:,1)*dum_pheno_std(cur_stream,:))  + grid_is(cur_stream,dum_int,pheno+1))
    end do
    dum_real=sum(dum_vect_grid(cur_stream,2:ninter_beta) + dum_vect_grid(cur_stream,1:(ninter_beta-1)))/(2*ninter_beta-2.)  !pour les bf: car prior uniforme (tester si pas uniforme)
    mean_isbf(mrk,pheno,1)=(mean_isbf(mrk,pheno,1)*real(iter-1)+ dum_real)/real(iter)
    mean_isbf(mrk,pheno,2)=(mean_isbf(mrk,pheno,2)*real(iter-1)+ dum_real**2)/real(iter)
    dum_vect_grid(cur_stream,1:(ninter_beta-1))=dum_vect_grid(cur_stream,2:ninter_beta) + dum_vect_grid(cur_stream,1:(ninter_beta-1))
    dum_real=(sum( dum_vect_grid(cur_stream,1:(ninter_beta-1)) * 0.5 * (grid_is(cur_stream,2:ninter_beta,1) + grid_is(cur_stream,1:(ninter_beta-1),1)) ))/sum(dum_vect_grid(cur_stream,:))
    mean_isbf(mrk,pheno,3)=(mean_isbf(mrk,pheno,3)*real(iter-1)+ dum_real)/real(iter)
    mean_isbf(mrk,pheno,4)=(mean_isbf(mrk,pheno,4)*real(iter-1)+ dum_real**2)/real(iter)
   end do
!$OMP END do   
!$OMP END PARALLEL 
  end do
 end if

end subroutine void_update_summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! IMPRIMER SUMMARY STATS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine void_print_summary()
!~  if(up_alpha_coop) then
!~   write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate Pij  = ',sum(acc_pij)/(real(thin)*real(nvaleurs)*real(nmrk)),&
!~                                           ' mean delta_pij= ',sum(delta_pij)/real(nmrk)
!~  else 
!~   write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate Pij  = ',sum(acc_pij)/(real(thin)*real(nvaleurs)*real(nmrk)*real(npop)),&
!~                                           ' mean delta_pij= ',sum(delta_pij)/(real(nmrk)*real(npop))
!~  end if                                        
!~  if(nich_prior) then
!~   write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate Cj   = ',sum(acc_nichc)/(real(thin)*real(nvaleurs)*real(npop)),&
!~                                         ' mean delta_cj= ',sum(delta_nichc)/real(npop)
!~  end if
!~  write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate Pi   = ',sum(acc_pi)/(real(thin)*real(nvaleurs)*real(nmrk)),&
!~                                         ' mean delta_pi= ',sum(delta_pi)/real(nmrk)
!~  if(up_params_beta) then 
!~   write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate b_mu = ',acc_beta_params(1)/(real(thin)*real(nvaleurs)),&
!~                                          '      delta  = ',delta_pi_beta_mu
!~   write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate b_phi= ',acc_beta_params(2)/(real(thin)*real(nvaleurs)),&
!~                                          '      delta  = ',delta_pi_beta_phi
!~  end if
 
 !!!!!!!!!!!
 !!impression summary stats
 !!!!!!!!!!!


 if(estim_omega) then
  open(3,file=sum_omega_file,status='unknown')
  if(nich_prior) then
   write (3,*) 'POP M_cj SD_cj DELTA_C ACC_C'
  else
   write (3,*) 'I J M_lda_ij SD_lda_ij M_omega_ij SD_omega_ij'
  end if
  open(4,file=sum_mat_omega_file,status='unknown') 
  do pop=1,npop 
   write(4,'(1000(f12.6,1x))') ( mean_omega(pop,dum_int,1), dum_int=1,npop )
   if(nich_prior) then
     write(3,'((i5,1x),4(f15.8,1x))') pop,mean_omega(pop,pop,1),sqrt(mean_omega(pop,pop,2)-(mean_omega(pop,pop,1))**2),&
									  delta_nichc(pop),acc_nichc(pop)/(real(nvaleurs)*real(thin))
   else
    do mrk=1,npop
     write(3,'(2(i5,1x),4(f15.8,1x))') pop,mrk,mean_lda(pop,mrk,1),sqrt(mean_lda(pop,mrk,2)-(mean_lda(pop,mrk,1))**2),&
                                               mean_omega(pop,mrk,1),sqrt(mean_omega(pop,mrk,2)-(mean_omega(pop,mrk,1))**2)
    end do
   end if
  end do
 end if
 close(3) ; close(4)

!Impression des données pop et locus
 if(poolseq) then
   open(1,file=sum_yij_pij_file,status='unknown')
!~    write (1,*) 'POP MRK M_Y SD_Y M_P SD_P M_Pstd SD_Pstd DELTA_P ACC_P DELTA_Y ACC_Y'
   write (1,*) 'POP MRK M_Y SD_Y M_P SD_P M_Pstd SD_Pstd'
 else
   open(1,file=sum_pij_file,status='unknown')
!~    write (1,*) 'POP MRK M_P SD_P M_Pstd SD_Pstd DELTA_P ACC_P'
   write (1,*) 'POP MRK M_P SD_P M_Pstd SD_Pstd'
 end if
 
 open(2,file=sum_pi_file,status='unknown') 
 if(root_pool) then
!~   write (2,*) 'MRK M_P SD_P DELTA_P ACC_P M_XtX SD_XtX XtXst log10(1/pval) M_Y SD_Y DELTA_Y ACC_Y' 
  write (2,*) 'MRK M_P SD_P M_XtX SD_XtX XtXst log10(1/pval) M_Y SD_Y' 
 else
!~   write (2,*) 'MRK M_P SD_P DELTA_P ACC_P M_XtX SD_XtX XtXst log10(1/pval)' 
  write (2,*) 'MRK M_P SD_P M_XtX SD_XtX XtXst log10(1/pval)' 
 end if
 dum_int=npop ; if(up_alpha_coop) dum_int=1
 !standardisation des XtX*: à parir des posterior Pstd
 dum_real=sum(mean_pij(:,:,3))/(real(nmrk)*real(npop))
 dum_real2=sqrt(sum( (mean_pij(:,:,3)-dum_real)**2)/(real(nmrk)*real(npop)-1.))
 dum_freq(:,:,1)=(mean_pij(:,:,3) - dum_real)/dum_real2
 do mrk=1,nmrk
  mean_xtx(mrk,3)=dot_product(dum_freq(mrk,:,1),dum_freq(mrk,:,1))
  call chi_square_cdf(mean_xtx(mrk,3),real(npop),dum_real)
  if(root_pool) then
   write(2,'(i8,1x,8(f12.8,1x))') mrk,mean_pi(mrk,1),sqrt(mean_pi(mrk,2)-(mean_pi(mrk,1))**2),&
!~    write(2,'(1i8,10(f12.8,1x),i5,1x,f6.4)') mrk,mean_pi(mrk,1),sqrt(mean_pi(mrk,2)-(mean_pi(mrk,1))**2),&
!~                                 delta_pi(mrk),acc_pi(mrk)/(real(nvaleurs)*real(thin)),&
                                mean_xtx(mrk,1),sqrt(mean_xtx(mrk,2)-(mean_xtx(mrk,1))**2),&
!~                                 mean_xtx(mrk,3),sqrt(mean_xtx(mrk,4)-(mean_xtx(mrk,3))**2),&
                                mean_xtx(mrk,3),-1*log10(1.-dum_real),&
                                mean_rooty(mrk,1),sqrt(mean_rooty(mrk,2)-(mean_rooty(mrk,1))**2)!,&
!~                                 delta_rooty(mrk),acc_rooty(mrk)/(real(nvaleurs)*real(thin)) 
  else
   write(2,'(i8,1x,6(f12.8,1x))') mrk,mean_pi(mrk,1),sqrt(mean_pi(mrk,2)-(mean_pi(mrk,1))**2),&
!~    write(2,'(1i8,8(f12.8,1x))') mrk,mean_pi(mrk,1),sqrt(mean_pi(mrk,2)-(mean_pi(mrk,1))**2),&
!~                                 delta_pi(mrk),acc_pi(mrk)/(real(nvaleurs)*real(thin)),&
                                mean_xtx(mrk,1),sqrt(mean_xtx(mrk,2)-(mean_xtx(mrk,1))**2),&
                                mean_xtx(mrk,3),-1*log10(1.-dum_real)
!                                mean_xtx(mrk,3),sqrt(mean_xtx(mrk,4)-(mean_xtx(mrk,3))**2)
  end if
  do pop=1,npop
   dum_int=pop ; if(up_alpha_coop) dum_int=1
   if(poolseq) then
     write(1,'(i4,1x,i8,1x,2(f8.2,1x),4(f12.8,1x))') pop,mrk,mean_yij(mrk,pop,1),sqrt(mean_yij(mrk,pop,2)-(mean_yij(mrk,pop,1))**2),&
!~      write(1,'(i4,1x,i8,2(f8.2,1x),6(f12.8,1x),1x,i5,1x,f6.4)') pop,mrk,mean_yij(mrk,pop,1),sqrt(mean_yij(mrk,pop,2)-(mean_yij(mrk,pop,1))**2),&
                                        mean_pij(mrk,pop,1),sqrt(mean_pij(mrk,pop,2)-(mean_pij(mrk,pop,1))**2),&
                                        mean_pij(mrk,pop,3),sqrt(mean_pij(mrk,pop,4)-(mean_pij(mrk,pop,3))**2)!,&
!~                                         delta_pij(mrk,dum_int),acc_pij(mrk,dum_int)/(real(nvaleurs)*real(thin)),&
!~                                         delta_y(mrk,pop),acc_y(mrk,pop)/(real(nvaleurs)*real(thin))
   else
     write(1,'(i4,1x,i8,1x,4(f12.8,1x))') pop,mrk,mean_pij(mrk,pop,1),sqrt(mean_pij(mrk,pop,2)-(mean_pij(mrk,pop,1))**2),&
!~      write(1,'(i4,1x,i8,6(f12.8,1x))') pop,mrk,mean_pij(mrk,pop,1),sqrt(mean_pij(mrk,pop,2)-(mean_pij(mrk,pop,1))**2),&
                                        mean_pij(mrk,pop,3),sqrt(mean_pij(mrk,pop,4)-(mean_pij(mrk,pop,3))**2)!,&
!~                                         delta_pij(mrk,dum_int),acc_pij(mrk,dum_int)/(real(nvaleurs)*real(thin))
   end if
  end do
 end do
  close(1) ; close(2)

  if(up_params_beta) then
    open(4,file=sum_betaparams_file,status='unknown') ;  write (4,*) 'PARAM Mean SD'
    write(4,'(A,1x,2(f12.6,1x))') 'a_beta_pi ',mean_beta_params(1,1),sqrt(mean_beta_params(1,2)-(mean_beta_params(1,1))**2)
    write(4,'(A,1x,2(f12.6,1x))') 'b_beta_pi ',mean_beta_params(2,1),sqrt(mean_beta_params(2,2)-(mean_beta_params(2,1))**2)
    close(4)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Impression des contrastes (le cas echeant)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(opt_contrast) then
  open(4,file=sum_contrast_file,status='unknown') 
!~   write (4,*) 'CONTRAST MRK M_C SD_C M_C2 SD_C2 C_std C_std2 M_Calt Calt_std2' !Calt calculé sur les freq vraies cad in (0,1)
  write (4,*) 'CONTRAST MRK M_C2 SD_C2 C2_std log10(1/pval)'! C2_std_old' !M_C2=C2 posterior; pval calulee sur le C2 standardise
  do contrast=1,ncontrast
!   dum_freq(:,1,2)=mean_contrast(contrast,:,1)
!   dum_real=sum(dum_freq(:,1,2))/real(nmrk)
!   dum_real2=sqrt(sum((dum_freq(:,1,2)-dum_real)**2)/real(nmrk-1))
!   dum_freq(:,1,2)=(dum_freq(:,1,2)-dum_real)/dum_real2
!~    dum_freq(:,1,2)=mean_contrast(contrast,:,5)
!~    dum_real=sum(dum_freq(:,1,2))/real(nmrk)
!~    dum_real2=sqrt(sum((dum_freq(:,1,2)-dum_real)**2)/real(nmrk-1))
!~    dum_freq(:,1,2)=(dum_freq(:,1,2)-dum_real)/dum_real2   
   do mrk=1,nmrk
    mean_contrast(contrast,mrk,5)=(sum(dum_freq(mrk,:,1)*contrastes(contrast,:)) / contrast_sd(contrast))**2
!~     write(4,'(i3,1x,i7,1x,8(f12.8,1x))') contrast,mrk,mean_contrast(contrast,mrk,1),sqrt(mean_contrast(contrast,mrk,2)-(mean_contrast(contrast,mrk,1))**2),&
!~                                           mean_contrast(contrast,mrk,3),sqrt(mean_contrast(contrast,mrk,4)-(mean_contrast(contrast,mrk,3))**2),&
!~                                           dum_freq(mrk,1,1),dum_freq(mrk,1,1)**2,mean_contrast(contrast,mrk,1),dum_freq(mrk,1,2)**2
    call chi_square_cdf(mean_contrast(contrast,mrk,5),1.0,dum_real)
    write(4,'(i3,1x,i7,1x,4(f12.8,1x))') contrast,mrk,mean_contrast(contrast,mrk,3),sqrt(mean_contrast(contrast,mrk,4)-(mean_contrast(contrast,mrk,3))**2),&
                                                  mean_contrast(contrast,mrk,5),-1*log10(1.-dum_real)!,dum_freq(mrk,1,2)**2
   end do
  end do
  close(4)
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!cas avec covariables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if(opt_pheno) then
  if(.not. estim_tau_beta) then 
   mean_tau_beta(:,1)=tau_beta0 ; mean_tau_beta(:,2)=tau_beta0**2 
  end if
  if(opt_aux) then !Bayes Factor
   open(25,file=sum_betai_file,status='unknown') !on ne met pas les accRate car pas de sens
   write (25,*) 'COVARIABLE MRK M_Beta SD_Beta PIP BF(dB) '
   open(26,file=sum_Pdelta_file,status='unknown')
   write (26,*) 'COVARIABLE M_P SD_P'!' M_Tau SD_tau'
   dum_real=P_beta_params(1)/sum(P_beta_params) ; dum_real=(1.-dum_real)/dum_real !inverse prior odd
!~    dum_real=(1.-mean_p_delta(pheno,1))/mean_p_delta(pheno,1) 
   do pheno=1,npheno
    write(26,'(i3,1x,4(f12.8,1x))') pheno,mean_p_delta(pheno,1),sqrt(mean_p_delta(pheno,2)-mean_p_delta(pheno,1)**2)!,mean_tau_beta(pheno,1),sqrt(mean_tau_beta(pheno,2)-mean_tau_beta(pheno,1)**2)
    do mrk=1,nmrk
     tmp_mean=mean_delta_i(mrk,pheno)/(1.-mean_delta_i(mrk,pheno)) !posterior odds
     if(mean_delta_i(mrk,pheno)<1./real(nvaleurs)) tmp_mean=0.5/(real(nvaleurs)-0.5)
     if(mean_delta_i(mrk,pheno)>(nvaleurs-1.)/nvaleurs) tmp_mean=2.*(real(nvaleurs)-0.5)
     bf=log10(dum_real*tmp_mean)
     write(25,'(i3,1x,i7,1x,10(f12.8,1x))') pheno,mrk,mean_beta_i(mrk,pheno,1),&
                                          sqrt(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2),&
                                          mean_delta_i(mrk,pheno),10.*bf
    end do
   end do
   close(26) ; close(25)
  else !KLD et BPval
   if(up_betai_coop) then
    open(25,file=sum_betai_file,status='unknown')
    write (25,*) 'COVARIABLE MRK M_Beta SD_Beta AccRateB DeltaB eBPmc'
   else
    open(25,file=sum_betai_file,status='unknown')
    write (25,*) 'COVARIABLE MRK M_Beta SD_Beta KLD logBPval'
   end if
   if(estim_tau_beta) then
    open(26,file=sum_taufile,status='unknown')
    write (26,*) 'COVARIABLE M_Tau SD_tau'
   end if
   do pheno=1,npheno
    if(estim_tau_beta) then
      write(26,'(i3,1x,4(f12.8,1x))') pheno,mean_tau_beta(pheno,1),sqrt(mean_tau_beta(pheno,2)-mean_tau_beta(pheno,1)**2)
    end if
    do mrk=1,nmrk
     bpval=mean_beta_i(mrk,pheno,1)/sqrt(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2)
     bpval=min(8.1,max(bpval,-8.1))
     call normal_01_cdf(bpval,dum_real)
     bpval=-1.*log10(1.-2.*abs(0.5-dum_real))
     if(up_betai_coop) then
      write(25,'(i3,1x,i7,1x,10(f12.8,1x))') pheno,mrk,mean_beta_i(mrk,pheno,1),&
                                           acc_betai(mrk,pheno)/(real(nvaleurs)*real(thin)),delta_betai(mrk,pheno),&
                                           sqrt(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2),bpval
     else
      dum_real=(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2)*mean_tau_beta(pheno,1)
      kld=-1.*log(dum_real) + dum_real
      dum_real= (mean_beta_i(mrk,pheno,2)**2)*mean_tau_beta(pheno,1)
      kld=kld+dum_real-1. ; kld=0.5*kld
      write(25,'(i3,1x,i7,1x,8(f12.8,1x))') pheno,mrk,mean_beta_i(mrk,pheno,1),&
                                           sqrt(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2),&
                                           kld,bpval
     end if
    end do
   end do
   close(25)
   if(estim_tau_beta) close(26)
  end if

  if(opt_predcov) then
   open(125,file=sum_covpred_file,status='unknown')
   write(125,*) 'COVARIABLE POP M_cov SD_cov'
   do pheno=1,npheno
    do pop=1,npop
      write(125,'(i3,1x,i4,1x,2(f12.8,1x))') pheno,pop,mean_phenoval(pop,pheno,1),sqrt(mean_phenoval(pop,pheno,2)-(mean_phenoval(pop,pheno,1))**2)
    end do
   end do
  end if
  
 end if

 if(opt_reg) then
  open(25,file=sum_beta_i_reg_file,status='unknown')
  write (25,*) 'COVARIABLE MRK M_Pearson SD_Pearson BF(dB) Beta_is SD_Beta_is eBPis' !pas de sens de mettre SD_BF
  do pheno=1,npheno
   do mrk=1,nmrk
    bpval=mean_isbf(mrk,pheno,3)/sqrt(mean_isbf(mrk,pheno,4)-(mean_isbf(mrk,pheno,3))**2)
    bpval=min(8.1,max(bpval,-8.1))
    call normal_01_cdf(bpval,dum_real)
    bpval=-1.*log10(1.-2.*abs(0.5-dum_real)) 
    write(25,'(i3,1x,i7,1x,8(f14.8,1x))') pheno,mrk,mean_rho_coef(mrk,pheno,1),&
                                       sqrt(mean_rho_coef(mrk,pheno,2)-(mean_rho_coef(mrk,pheno,1))**2),&
                                       10.*log10(mean_isbf(mrk,pheno,1)),&
                                       mean_isbf(mrk,pheno,3),sqrt(mean_isbf(mrk,pheno,4)-(mean_isbf(mrk,pheno,3))**2),bpval
   end do
  end do
  close(25)
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !impression deviance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(.not. indseq) then
 lpml=0.
 do mrk=1,nmrk
  do pop=1,npop
   lpml=lpml - log(CPO(mrk,pop))
  end do
 end do
 open(1002,file=sum_DIC_file,status='unknown')
 write (1002,*) 'bar(D) pD DIC LPML'
 dum_real16=deviance(Y_OBS,N_OBS,mean_pij(:,:,1),log_cn)
 dum_real16_2=mean_deviance-dum_real16
 write(1002,'(4(f20.2,1x))') mean_deviance,dum_real16_2,dum_real16+2.*dum_real16_2,lpml
 close(1002)
end if

end subroutine void_print_summary

end program baypass



