module covid19abm
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
include("matrices_code.jl")
@enum HEALTH SUS LAT PRE ASYMP MILD MISO INF IISO HOS ICU REC DED  LAT2 PRE2 ASYMP2 MILD2 MISO2 INF2 IISO2 HOS2 ICU2 REC2 DED2 LAT3 PRE3 ASYMP3 MILD3 MISO3 INF3 IISO3 HOS3 ICU3 REC3 DED3 UNDEF 

Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    swap::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    wentTo::HEALTH = UNDEF
    sickby::Int64 = -1
    nextday_meetcnt::Int16 = 0 ## how many contacts for a single day
    age::Int16   = 0    # in years. don't really need this but left it incase needed later
    ag::Int16   = 0
    tis::Int16   = 0   # time in state 
    exp::Int16   = 0   # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16   = 999   # day of infection.
    iso::Bool = false  ## isolated (limited contacts)
    isovia::Symbol = :null ## isolated via quarantine (:qu), preiso (:pi), intervention measure (:im), or contact tracing (:ct)    
    tracing::Bool = false ## are we tracing contacts for this individual?
    tracestart::Int16 = -1 ## when to start tracing, based on values sampled for x.dur
    traceend::Int16 = -1 ## when to end tracing
    tracedby::UInt32 = 0 ## is the individual traced? property represents the index of the infectious person 
    tracedxp::Int16 = 0 ## the trace is killed after tracedxp amount of days
    comorbidity::Int8 = 0 ##does the individual has any comorbidity?
    vac_status::Int8 = 0 ##
    vac_ef_symp::Float16 = 0.0 
    vac_ef_inf::Float16 = 0.0 
    vac_ef_sev::Float16 = 0.0

    got_inf::Bool = false
    herd_im::Bool = false
    hospicu::Int8 = -1
    ag_new::Int16 = -1
    hcw::Bool = false
    days_vac::Int64 = -1
    vac_red::Float64 = 0.0
    first_one::Bool = false
    strain::Int16 = -1
    index_day::Int64 = 1
    relaxed::Bool = false
    recovered::Bool = false
    vaccine::Symbol = :none
    vaccine_n::Int16 = 0
    protected::Int64 = 0
end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0345       
    seasonal::Bool = false ## seasonal betas or not
    popsize::Int64 = 100000
    prov::Symbol = :usa
    calibration::Bool = false
    calibration2::Bool = false 
    heatmap::Bool = false
    ignore_cal::Bool = false
    start_several_inf::Bool = true
    modeltime::Int64 = 500
    initialinf::Int64 = 50
    initialhi::Int64 = 0 ## initial herd immunity, inserts number of REC individuals
    τmild::Int64 = 0 ## days before they self-isolate for mild cases
    fmild::Float64 = 0.0  ## percent of people practice self-isolation
    fsevere::Float64 = 0.0 #
    eldq::Float64 = 0.0 ## complete isolation of elderly
    eldqag::Int8 = 5 ## default age group, if quarantined(isolated) is ag 5. 
    fpreiso::Float64 = 0.0 ## percent that is isolated at the presymptomatic stage
    tpreiso::Int64 = 0## preiso is only turned on at this time. 
    frelasymp::Float64 = 0.26 ## relative transmission of asymptomatic
    ctstrat::Int8 = 0 ## strategy 
    fctcapture::Float16 = 0.0 ## how many symptomatic people identified
    fcontactst::Float16 = 0.0 ## fraction of contacts being isolated/quarantined
    cidtime::Int8 = 0  ## time to identification (for CT) post symptom onset
    cdaysback::Int8 = 0 ## number of days to go back and collect contacts
    #vaccine_ef::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    vac_com_dec_max::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    vac_com_dec_min::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    herd::Int8 = 0 #typemax(Int32) ~ millions
    set_g_cov::Bool = false ###Given proportion for coverage
    cov_val::Float64 = 0.2
    file_index::Int16 = 0
    
    #the cap for coverage should be 90% for 65+; 95% for HCW; 80% for 50-64; 60% for 16-49; and then 50% for 12-15 (starting from June 1).

    hcw_vac_comp::Float64 = 0.95
    hcw_prop::Float64 = 0.05 #prop que é trabalhador da saude
    
    eld_comp::Float64 = 0.95
    old_adults::Float64 = 0.95
    young_adults::Float64 = 0.8
    kid_comp::Float64 = 0.8
    #comor_comp::Float64 = 0.7 #prop comorbidade tomam

    vac_period::Array{Int64,1} = [21;28]
    n_comor_comp::Float64 = 1.0
    min_age_vac::Int64 = 18
    
    #=------------ Vaccine Efficacy ----------------------------=#
    days_to_protection::Array{Array{Int64,1},1} = [[14],[0;7]]
    vac_efficacy_inf::Array{Array{Array{Float64,1},1},1} = [[[0.46],[0.6;0.92]],[[0.295],[0.6;0.895]]] #### 50:5:80
    vac_efficacy_symp::Array{Array{Array{Float64,1},1},1} = [[[0.57],[0.66;0.94]],[[0.536],[0.62;0.934]]] #### 50:5:80
    vac_efficacy_sev::Array{Array{Array{Float64,1},1},1} = [[[0.62],[0.80;0.92]],[[0.541],[0.8;0.94]]]#### 50:5:80


    vaccinating::Bool = false #vaccinating?
    single_dose::Bool = false #unique dose
    drop_rate::Float64 = 0.0 #probability of not getting second dose
    fixed_cov::Float64 = 0.4 #coverage of population
    day_insert_kids::Int64 = 255


    pfizer_proportion::Float64 = 1.0
    red_risk_perc::Float64 = 1.0 #relative isolation in vaccinated individuals
    reduction_protection::Float64 = 0.0 #reduction in protection against infection
    extra_dose::Bool = false #if it is turned on, extradoses are given to general population
    extra_dose_n::Int64 = 0 #how many
    extra_dose_day::Int64 = 999 #when extra doses are implemented
    days_Rt::Array{Int64,1} = [100;200;300] #days to get Rt

    sec_strain_trans::Float64 = 1.5#1.5 #transmissibility of second strain
    ins_sec_strain::Bool = true #insert second strain?
    initialinf2::Int64 = 1 #number of initial infected of second strain
    time_sec_strain::Int64 = 92 #when will the second strain introduced

    ins_third_strain::Bool = false #insert third strain?
    initialinf3::Int64 = 5 #number of initial infected of third strain
    time_third_strain::Int64 = 150 #when will the third strain introduced
    third_strain_trans::Float64 = 1.0 #transmissibility of third strain
    reduction_recovered::Float64 = 0.21

    max_vac_delay::Int64 = 42 #max delay before protection starts waning
    min_eff = 0.02 #min efficacy when waning
    vac_effect::Int64 = 1 #vac effect, if 1 the difference between doses is added to first, if 2 the second dose is always vac_efficacy
    no_cap::Bool = true ## no maximum coverage
    mortality_inc::Float64 = 1.3 #The mortality increase when infected by strain 2

    time_change::Int64 = 999## used to calibrate the model
    how_long::Int64 = 1## used to calibrate the model
    how_much::Float64 = 0.0## used to calibrate the model
    rate_increase::Float64 = how_much/how_long## used to calibrate the model
    time_change_contact::Array{Int64,1} = [1;map(y->47+y,0:(5));map(y->81+y,0:(5));map(y->93+y,0:4);map(y->111+y,0:8);map(y->127+y,0:11);map(y->180+y,0:(9))]
    change_rate_values::Array{Float64,1} = [1;map(y->1.0+0.01*y,1:6);map(y->1.06-0.01*y,1:6);map(y->1.0-(0.03/5)*y,1:5);map(y->0.97+(0.084/9)*y,1:9);map(y->1.054-(0.20/12)*y,1:12);map(y->0.864+(0.045/10)*y,1:10)]
    contact_change_rate::Float64 = 1.0 #the rate that receives the value of change_rate_values
    contact_change_2::Float64 = 0.5 ##baseline number that multiplies the contact rate

    relaxed::Bool = false
    relaxing_time::Int64 = 999#215 ### relax measures for vaccinated
    status_relax::Int16 = 2
    relax_after::Int64 = 1

    priority::Bool = true
    scenario::Symbol = :statuscuo

    time_back_to_normal::Int64 = 999 ###relaxing time of measures for non-vaccinated
    ### after calibration, how much do we want to increase the contact rate... in this case, to reach 70%
    ### 0.5*0.95 = 0.475, so we want to multiply this by 1.473684211
    back_normal_rate::Float64 = 1.506996771 #2.15285253=>1 # ####1.506996771 =>0.7 ###1.82992465 =>0.85
    
end

Base.@kwdef mutable struct ct_data_collect
    total_symp_id::Int64 = 0  # total symptomatic identified
    totaltrace::Int64 = 0     # total contacts traced
    totalisolated::Int64 = 0  # total number of people isolated
    iso_sus::Int64 = 0        # total susceptible isolated 
    iso_lat::Int64 = 0        # total latent isolated
    iso_asymp::Int64 = 0      # total asymp isolated
    iso_symp::Int64 = 0       # total symp (mild, inf) isolated
end

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## constants 
const humans = Array{Human}(undef, 0) 
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]
const agebraks_vac = @SVector [12:15, 16:17, 18:24, 25:39, 40:49, 50:64, 65:74, 75:99]
const BETAS = Array{Float64, 1}(undef, 0) ## to hold betas (whether fixed or seasonal), array will get resized
const ct_data = ct_data_collect()
export ModelParameters, HEALTH, Human, humans, BETAS

function runsim(simnum, ip::ModelParameters)
    # function runs the `main` function, and collects the data as dataframes. 
    hmatrix,hh1,hh2,hh3 = main(ip,simnum)            
    # get infectors counters
    infectors = _count_infectors()

    ct_numbers = (ct_data.total_symp_id, ct_data.totaltrace, ct_data.totalisolated, 
                    ct_data.iso_sus, ct_data.iso_lat, ct_data.iso_asymp, ct_data.iso_symp)
    ###use here to create the vector of comorbidity
    # get simulation age groups
    #ags = [x.ag for x in humans] # store a vector of the age group distribution 
    ags = [x.ag_new for x in humans] # store a vector of the age group distribution 
    all = _collectdf(hmatrix)
    spl = _splitstate(hmatrix, ags)
    ag1 = _collectdf(spl[1])
    ag2 = _collectdf(spl[2])
    ag3 = _collectdf(spl[3])
    ag4 = _collectdf(spl[4])
    ag5 = _collectdf(spl[5])
    ag6 = _collectdf(spl[6])
    insertcols!(all, 1, :sim => simnum); insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum);
    insertcols!(ag6, 1, :sim => simnum);

    
    R01 = zeros(Float64,size(hh1,1))

    for i = 1:size(hh1,1)
        if length(hh1[i]) > 0
            R01[i] = length(findall(k -> k.sickby in hh1[i],humans))/length(hh1[i])
        end
    end

    R02 = zeros(Float64,size(hh2,1))

    for i = 1:size(hh2,1)
        if length(hh2[i]) > 0
            R02[i] = length(findall(k -> k.sickby in hh2[i],humans))/length(hh2[i])
        end
    end

    R03 = zeros(Float64,size(hh3,1))

    for i = 1:size(hh3,1)
        if length(hh3[i]) > 0
            R03[i] = length(findall(k -> k.sickby in hh3[i],humans))/length(hh3[i])
        end
    end

    #return (a=all, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5, infectors=infectors, vi = vac_idx,ve=vac_ef_i,com = comorb_idx,n_vac = n_vac,n_inf_vac = n_inf_vac,n_inf_nvac = n_inf_nvac)
    return (a=all, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6, infectors=infectors,   
    iniiso = ct_data.totalisolated,
    R01 = R01,
    R02 = R02)
end
export runsim

function main(ip::ModelParameters,sim::Int64)
    Random.seed!(sim*726)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)  #logic: outside "ip" parameters are copied to internal "p" which is a global const and available everywhere. 

    p.popsize == 0 && error("no population size given")
    
    hmatrix = zeros(Int16, p.popsize, p.modeltime)
    initialize() # initialize population
    
     #h_init::Int64 = 0
    # insert initial infected agents into the model
    # and setup the right swap function. 
    if p.start_several_inf
        N = herd_immu_dist_4(sim,1)
        if p.initialinf > 0
            insert_infected(PRE, p.initialinf, 4, 1)[1]
        end
        #findall(x->x.health in (MILD,INF,LAT,PRE,ASYMP),humans)
    else
        #applying_vac(sim)
        herd_immu_dist_4(sim,1)
        insert_infected(PRE, 1, 4,1)[1]
    end    
    h_init1 = findall(x->x.health  in (LAT,MILD,MISO,INF,PRE,ASYMP),humans)
    h_init1 = [h_init1]
    h_init2 = []
    h_init3 = []
    ## save the preisolation isolation parameters
    
    p.fpreiso = 0

    # split population in agegroups 
    grps = get_ag_dist()
    count_change::Int64 = 1
   
    time_vac::Int64 = 1
    time_pos::Int64 = 0
    if p.vaccinating
        vac_ind = vac_selection(sim)
        
    else
        time_vac = 9999 #this guarantees that no one will be vaccinated
    end
    # start the time loop
    for st = 1:p.modeltime
        if p.ins_sec_strain && st == p.time_sec_strain ##insert second strain
            insert_infected(PRE, p.initialinf2, 4, 2)[1]
            h_init2 = findall(x->x.health  in (LAT2,MILD2,INF2,PRE2,ASYMP2),humans)
            h_init2 = [h_init2]
        end
        if p.ins_third_strain && st == p.time_third_strain #insert third strain
            insert_infected(PRE, p.initialinf3, 4, 3)[1]
            h_init3 = findall(x->x.health  in (LAT3,MILD3,INF3,PRE3,ASYMP3),humans)
            h_init3 = [h_init3]
        end
        if length(p.time_change_contact) >= count_change && p.time_change_contact[count_change] == st ###change contact pattern throughout the time
            setfield!(p, :contact_change_rate, p.change_rate_values[count_change])
            count_change += 1
        end
        # start of day
       
        if st == p.relaxing_time ### time that people vaccinated people is allowed to go back to normal
            setfield!(p, :relaxed, true)
        elseif st == p.time_back_to_normal ##time that non-vaccinated people is allowed to go back to normal
            setfield!(p, :contact_change_2, p.contact_change_2*p.back_normal_rate)
            #setfield!(p, :contact_change_rate, 1.0)
        end

        if time_pos < length(vaccination_days) && time_vac == vaccination_days[time_pos+1]
            time_pos += 1
        end
        time_vac += 1
        time_pos > 0 && vac_time!(sim,vac_ind,time_pos+1)
       
       
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps,sim)
        if st in p.days_Rt ### saves individuals that became latent on days_Rt
            aux1 = findall(x->x.swap == LAT,humans)
            h_init1 = vcat(h_init1,[aux1])
            aux2 = findall(x->x.swap == LAT2,humans)
            h_init2 = vcat(h_init2,[aux2])
            aux3 = findall(x->x.swap == LAT3,humans)
            h_init3 = vcat(h_init3,[aux3])
        end
        sw = time_update() ###update the system
        # end of day
    end
    
    
    return hmatrix,h_init1,h_init2,h_init3 ## return the model state as well as the age groups. 
end
export main


function vac_selection(sim::Int64)
    

    if p.priority
        aux_1 = map(k-> findall(y-> y.age in k && y.comorbidity == 1,humans),agebraks_vac)
        aux_2 = map(k-> findall(y-> y.age in k && y.comorbidity == 0,humans),agebraks_vac)

        v = map(x-> [aux_1[x];aux_2[x]],1:length(aux_1))
    else
        v = map(k-> findall(y-> y.age in k,humans),agebraks_vac)

    end
    

end


function vac_time!(sim::Int64,vac_ind::Vector{Vector{Int64}},time_pos::Int64)
    aux_states = (MILD, MISO, INF, IISO, HOS, ICU, DED,MILD2, MISO2, INF2, IISO2, HOS2, ICU2, DED2,MILD3, MISO3, INF3, IISO3, HOS3, ICU3, DED3)
    ##first dose
    rng = MersenneTwister(123*sim)
    ### lets create distribute the number of doses per age group
    
    remaining_doses::Int64 = 0

    for i in 1:length(vac_ind)
        pos = findall(y-> humans[y].vac_status == 1 && humans[y].days_vac >= p.vac_period[humans[y].vaccine_n] && !(humans[y].health in aux_states),vac_ind[i])
        
        l1 = min(vac_rate_2[time_pos,i],length(pos))
        remaining_doses += vac_rate_2[time_pos,i] - l1
        for j = 1:l1
            x = humans[vac_ind[i][pos[j]]]
            x.days_vac = 0
            x.vac_status = 2
            x.index_day = 1
        end

        pos = findall(y-> humans[y].vac_status == 0 && !(humans[y].health in aux_states),vac_ind[i])
        
        l2 = min(vac_rate_1[time_pos,i],length(pos))

        remaining_doses += vac_rate_1[time_pos,i] - l2

        for j = 1:l2
            x = humans[vac_ind[i][pos[j]]]
            x.days_vac = 0
            x.vac_status = 1
            x.index_day = 1
    
            x.vaccine = rand(rng) <= p.pfizer_proportion ? :pfizer : :moderna
            x.vaccine_n = x.vaccine == :pfizer ? 1 : 2
        end

    end

    ###remaining_doses are given to any individual within the groups that are vaccinated on that day
    if remaining_doses > 0

        pos = map(k->findall(y-> humans[y].vac_status == 1 && humans[y].days_vac >= p.vac_period[humans[y].vaccine_n] && !(humans[y].health in aux_states),vac_ind[k]),1:length(vac_ind))
        pos2 = map(k->findall(y-> humans[y].vac_status == 0 && !(humans[y].health in aux_states),vac_ind[k]),1:length(vac_ind))
        
        aux = findall(x-> vac_rate_1[time_pos,x] > 0 || vac_rate_2[time_pos,x] > 0, 1:length(vac_ind))
        position = map(k-> vac_ind[k][pos[k]],aux)
        position2 = map(k-> vac_ind[k][pos2[k]],aux)

        r = vcat(position...,position2...)
        m = min(remaining_doses,length(r))

        rr = sample(rng,r,m,replace=false)
        
        for i in rr
            x = humans[i]
            if x.vac_status == 0
                x.days_vac = 0
                x.vac_status = 1
                x.index_day = 1
                x.vaccine = rand(rng) <= p.pfizer_proportion ? :pfizer : :moderna
                x.vaccine_n = x.vaccine == :pfizer ? 1 : 2
            elseif x.vac_status == 1
                x.days_vac = 0
                x.vac_status = 2
                x.index_day = 1
            else
                error("error in humans vac status - vac time")
            end
        end
    end

end

function vac_update(x::Human)
   
    if x.vac_status == 1
        #x.index_day == 2 && error("saiu com indice 2")
        if x.days_vac == p.days_to_protection[x.vac_status][1]#14
            x.protected = 1
            x.index_day = min(length(p.days_to_protection[x.vac_status]),x.index_day+1)
        elseif x.days_vac == p.days_to_protection[x.vac_status][x.index_day]#14
            x.protected = x.index_day
            x.index_day = min(length(p.days_to_protection[x.vac_status]),x.index_day+1)
        end
        if !x.relaxed
            x.relaxed = p.relaxed &&  x.vac_status >= p.status_relax && x.days_vac >= p.relax_after ? true : false
        end
        x.days_vac += 1

    elseif x.vac_status == 2
        if x.days_vac == p.days_to_protection[x.vac_status][1]#0
            x.protected = 1
            x.index_day = min(length(p.days_to_protection[x.vac_status]),x.index_day+1)

        elseif x.days_vac == p.days_to_protection[x.vac_status][x.index_day]#7
            x.protected = x.index_day
            x.index_day = min(length(p.days_to_protection[x.vac_status]),x.index_day+1)
        end
        if !x.relaxed
            x.relaxed = p.relaxed &&  x.vac_status >= p.status_relax && x.days_vac >= p.relax_after ? true : false
        end
        x.days_vac += 1
    end
   
end
function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end

    # reset the contact tracing data collection structure
    for x in propertynames(ct_data)
        setfield!(ct_data, x, 0)
    end

    # resize and update the BETAS constant array
    init_betas()

    # resize the human array to change population size
    resize!(humans, p.popsize)
end
export reset_params, reset_params_default

function _model_check() 
    ## checks model parameters before running 
    (p.fctcapture > 0 && p.fpreiso > 0) && error("Can not do contact tracing and ID/ISO of pre at the same time.")
    (p.fctcapture > 0 && p.maxtracedays == 0) && error("maxtracedays can not be zero")
end

## Data Collection/ Model State functions
function _get_model_state(st, hmatrix)
    # collects the model state (i.e. agent status at time st)
    for i=1:length(humans)
        hmatrix[i, st] = Int(humans[i].health)
    end    
end
export _get_model_state

function _collectdf(hmatrix)
    ## takes the output of the humans x time matrix and processes it into a dataframe
    #_names_inci = Symbol.(["lat_inc", "mild_inc", "miso_inc", "inf_inc", "iiso_inc", "hos_inc", "icu_inc", "rec_inc", "ded_inc"])    
    #_names_prev = Symbol.(["sus", "lat", "mild", "miso", "inf", "iiso", "hos", "icu", "rec", "ded"])
    mdf_inc, mdf_prev = _get_incidence_and_prev(hmatrix)
    mdf = hcat(mdf_inc, mdf_prev)    
    _names_inc = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_INC"))
    _names_prev = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_PREV"))
    _names = vcat(_names_inc..., _names_prev...)
    datf = DataFrame(mdf, _names)
    insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    return datf
end

function _splitstate(hmatrix, ags)
    #split the full hmatrix into 4 age groups based on ags (the array of age group of each agent)
    #sizes = [length(findall(x -> x == i, ags)) for i = 1:4]
    matx = []#Array{Array{Int64, 2}, 1}(undef, 4)
    for i = 1:maximum(ags)#length(agebraks)
        idx = findall(x -> x == i, ags)
        push!(matx, view(hmatrix, idx, :))
    end
    return matx
end
export _splitstate

function _get_incidence_and_prev(hmatrix)
    cols = instances(HEALTH)[1:end - 1] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(hmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(hmatrix, cols[i])
    end
    return inc, pre
end

function _get_column_incidence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for r in eachrow(hmatrix)
        idx = findfirst(x -> x == inth, r)
        if idx !== nothing 
            timevec[idx] += 1
        end
    end
    return timevec
end

function herd_immu_dist_4(sim::Int64,strain::Int64)
    rng = MersenneTwister(200*sim)
    vec_n = zeros(Int32,6)
    N::Int64 = 0
    if p.herd == 5
        vec_n = [9; 148; 262;  68; 4; 9]
        N = 5

    elseif p.herd == 10
        vec_n = [32; 279; 489; 143; 24; 33]

        N = 9

    elseif p.herd == 20
        vec_n = [71; 531; 962; 302; 57; 77]

        N = 14
    elseif p.herd == 30
        vec_n = [105; 757; 1448; 481; 87; 122]

        N = 16
    elseif p.herd == 50
        vec_n = map(y->y*5,[32; 279; 489; 143; 24; 33])

        N = 16
    elseif p.herd == 0
        vec_n = [0;0;0;0;0;0]
       
    else
        error("No herd immunity")
    end

    for g = 1:6
        pos = findall(y->y.ag_new == g && y.health == SUS,humans)
        n_dist = min(length(pos),Int(floor(vec_n[g]*p.popsize/10000)))
        pos2 = sample(rng,pos,n_dist,replace=false)
        for i = pos2
            humans[i].strain = strain
            humans[i].swap = strain == 1 ? REC : REC2
            move_to_recovered(humans[i])
            humans[i].sickfrom = INF
            humans[i].herd_im = true
        end
    end
    return N
end

function _get_column_prevalence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for (i, c) in enumerate(eachcol(hmatrix))
        idx = findall(x -> x == inth, c)
        if idx !== nothing
            ps = length(c[idx])    
            timevec[i] = ps    
        end
    end
    return timevec
end

function _count_infectors()     
    pre_ctr = asymp_ctr = mild_ctr = inf_ctr = pre_ctr2 = asymp_ctr2 = mild_ctr2 = inf_ctr2 = pre_ctr3 = asymp_ctr3 = mild_ctr3 = inf_ctr3 = 0
    for x in humans 
        if x.health != SUS ## meaning they got sick at some point
            if x.sickfrom == PRE
                pre_ctr += 1
            elseif x.sickfrom == ASYMP
                asymp_ctr += 1
            elseif x.sickfrom == MILD || x.sickfrom == MISO 
                mild_ctr += 1 
            elseif x.sickfrom == INF || x.sickfrom == IISO 
                inf_ctr += 1 
            elseif x.sickfrom == PRE2
                pre_ctr2 += 1
            elseif x.sickfrom == ASYMP2
                asymp_ctr2 += 1
            elseif x.sickfrom == MILD2 || x.sickfrom == MISO2 
                mild_ctr2 += 1 
            elseif x.sickfrom == INF2 || x.sickfrom == IISO2 
                inf_ctr2 += 1 
            elseif x.sickfrom == PRE3
                pre_ctr3 += 1
            elseif x.sickfrom == ASYMP3
                asymp_ctr3 += 1
            elseif x.sickfrom == MILD3 || x.sickfrom == MISO3 
                mild_ctr3 += 1 
            elseif x.sickfrom == INF3 || x.sickfrom == IISO3 
                inf_ctr3 += 1 
            else 
                error("sickfrom not set right: $(x.sickfrom)")
            end
        end
    end
    return (pre_ctr, asymp_ctr, mild_ctr, inf_ctr,pre_ctr2, asymp_ctr2, mild_ctr2, inf_ctr2,pre_ctr3, asymp_ctr3, mild_ctr3, inf_ctr3)
end

export _collectdf, _get_incidence_and_prev, _get_column_incidence, _get_column_prevalence, _count_infectors

## initialization functions 
function get_province_ag(prov) 
    ret = @match prov begin        
        #=:alberta => Distributions.Categorical(@SVector [0.0655, 0.1851, 0.4331, 0.1933, 0.1230])
        :bc => Distributions.Categorical(@SVector [0.0475, 0.1570, 0.3905, 0.2223, 0.1827])
        :canada => Distributions.Categorical(@SVector [0.0540, 0.1697, 0.3915, 0.2159, 0.1689])
        :manitoba => Distributions.Categorical(@SVector [0.0634, 0.1918, 0.3899, 0.1993, 0.1556])
        :newbruns => Distributions.Categorical(@SVector [0.0460, 0.1563, 0.3565, 0.2421, 0.1991])
        :newfdland => Distributions.Categorical(@SVector [0.0430, 0.1526, 0.3642, 0.2458, 0.1944])
        :nwterrito => Distributions.Categorical(@SVector [0.0747, 0.2026, 0.4511, 0.1946, 0.0770])
        :novasco => Distributions.Categorical(@SVector [0.0455, 0.1549, 0.3601, 0.2405, 0.1990])
        :nunavut => Distributions.Categorical(@SVector [0.1157, 0.2968, 0.4321, 0.1174, 0.0380])
        
        :pei => Distributions.Categorical(@SVector [0.0490, 0.1702, 0.3540, 0.2329, 0.1939])
        :quebec => Distributions.Categorical(@SVector [0.0545, 0.1615, 0.3782, 0.2227, 0.1831])
        :saskat => Distributions.Categorical(@SVector [0.0666, 0.1914, 0.3871, 0.1997, 0.1552])
        :yukon => Distributions.Categorical(@SVector [0.0597, 0.1694, 0.4179, 0.2343, 0.1187])=#
        :ontario => Distributions.Categorical(@SVector [0.0519, 0.1727, 0.3930, 0.2150, 0.1674])
        :usa => Distributions.Categorical(@SVector [0.059444636404977,0.188450296592341,0.396101793107413,0.189694011721906,0.166309262173363])
        :newyork   => Distributions.Categorical(@SVector [0.064000, 0.163000, 0.448000, 0.181000, 0.144000])
        _ => error("shame for not knowing your canadian provinces and territories")
    end       
    return ret  
end
export get_province_ag

function comorbidity(ag::Int16)

    a = [4;19;49;64;79;999]
    g = findfirst(x->x>=ag,a)
    prob = [0.05; 0.1; 0.28; 0.55; 0.74; 0.81]

    com = rand() < prob[g] ? 1 : 0

    return com    
end
export comorbidity


function initialize() 
    agedist = get_province_ag(p.prov)
    for i = 1:p.popsize 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        x.ag = rand(agedist)
        x.age = rand(agebraks[x.ag]) 
        a = [4;19;49;64;79;999]
        g = findfirst(y->y>=x.age,a)
        x.ag_new = g
        x.exp = 999  ## susceptible people don't expire.
        x.dur = sample_epi_durations() # sample epi periods   
        if rand() < p.eldq && x.ag == p.eldqag   ## check if elderly need to be quarantined.
            x.iso = true   
            x.isovia = :qu         
        end
        x.comorbidity = comorbidity(x.age)
        x.vac_red = p.vac_com_dec_min+rand()*(p.vac_com_dec_max-p.vac_com_dec_min)
        # initialize the next day counts (this is important in initialization since dyntrans runs first)
        get_nextday_counts(x)
        
    end
end
export initialize

function init_betas() 
    if p.seasonal  
        tmp = p.β .* td_seasonality()
    else 
        tmp = p.β .* ones(Float64, p.modeltime)
    end
    resize!(BETAS, length(tmp))
    for i = 1:length(tmp)
        BETAS[i] = tmp[i]
    end
end

function td_seasonality()
    ## returns a vector of seasonal oscillations
    t = 1:p.modeltime
    a0 = 6.261
    a1 = -11.81
    b1 = 1.817
    w = 0.022 #0.01815    
    temp = @. a0 + a1*cos((80-t)*w) + b1*sin((80-t)*w)  #100
    #temp = @. a0 + a1*cos((80-t+150)*w) + b1*sin((80-t+150)*w)  #100
    temp = (temp .- 2.5*minimum(temp))./(maximum(temp) .- minimum(temp)); # normalize  @2
    return temp
end

function get_ag_dist() 
    # splits the initialized human pop into its age groups
    grps =  map(x -> findall(y -> y.ag == x, humans), 1:length(agebraks)) 
    return grps
end

function insert_infected(health, num, ag,strain) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    l = findall(x -> x.health == SUS && x.ag == ag, humans)
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            x = humans[i]
            x.strain = strain
            x.first_one = true

            if x.strain == 1
                if health == PRE
                    x.swap =  PRE 
                    move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
                elseif health == LAT
                    x.swap = LAT 
                    move_to_latent(x)
                elseif health == MILD
                    x.swap =  MILD 
                    move_to_mild(x)
                elseif health == INF
                    x.swap = INF
                    move_to_infsimple(x)
                elseif health == ASYMP
                    x.swap = ASYMP 
                    move_to_asymp(x)
                elseif health == REC 
                    x.swap = REC
                    move_to_recovered(x)
                else 
                    error("can not insert human of health $(health)")
                end
            elseif x.strain == 2
                if health == PRE
                    x.swap =  PRE2 
                    move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
                elseif health == LAT
                    x.swap = LAT2 
                    move_to_latent(x)
                elseif health == MILD
                    x.swap =  MILD2 
                    move_to_mild(x)
                elseif health == INF
                    x.swap = INF2
                    move_to_infsimple(x)
                elseif health == ASYMP
                    x.swap = ASYMP2 
                    move_to_asymp(x)
                elseif health == REC 
                    x.swap = REC2
                    move_to_recovered(x)
                else 
                    error("can not insert human of health $(health)")
                end
            elseif x.strain == 3
                if health == PRE
                    x.swap =  PRE3 
                    move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
                elseif health == LAT
                    x.swap = LAT3 
                    move_to_latent(x)
                elseif health == MILD
                    x.swap =  MILD3 
                    move_to_mild(x)
                elseif health == INF
                    x.swap = INF3
                    move_to_infsimple(x)
                elseif health == ASYMP
                    x.swap = ASYMP3 
                    move_to_asymp(x)
                elseif health == REC 
                    x.swap = REC3
                    move_to_recovered(x)
                else 
                    error("can not insert human of health $(health)")
                end
            else
                error("no strain")
            end
            x.sickfrom = INF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
            
        end
    end    
    return h
end
export insert_infected

function time_update()
    # counters to calculate incidence
    lat=0; pre=0; asymp=0; mild=0; miso=0; inf=0; infiso=0; hos=0; icu=0; rec=0; ded=0;
    lat2=0; pre2=0; asymp2=0; mild2=0; miso2=0; inf2=0; infiso2=0; hos2=0; icu2=0; rec2=0; ded2=0;
    lat3=0; pre3=0; asymp3=0; mild3=0; miso3=0; inf3=0; infiso3=0; hos3=0; icu3=0; rec3=0; ded3=0;
    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
        if x.tis >= x.exp             
            @match Symbol(x.swap) begin
                :LAT  => begin move_to_latent(x); lat += 1; end
                :PRE  => begin move_to_pre(x); pre += 1; end
                :ASYMP => begin move_to_asymp(x); asymp += 1; end
                :MILD => begin move_to_mild(x); mild += 1; end
                :MISO => begin move_to_miso(x); miso += 1; end
                :INF  => begin move_to_inf(x); inf +=1; end    
                :IISO => begin move_to_iiso(x); infiso += 1; end
                :HOS  => begin move_to_hospicu(x); hos += 1; end 
                :ICU  => begin move_to_hospicu(x); icu += 1; end
                :REC  => begin move_to_recovered(x); rec += 1; end
                :DED  => begin move_to_dead(x); ded += 1; end
                :LAT2  => begin move_to_latent(x); lat2 += 1; end
                :PRE2  => begin move_to_pre(x); pre2 += 1; end
                :ASYMP2 => begin move_to_asymp(x); asymp2 += 1; end
                :MILD2 => begin move_to_mild(x); mild2 += 1; end
                :MISO2 => begin move_to_miso(x); miso2 += 1; end
                :INF2  => begin move_to_inf(x); inf2 +=1; end    
                :IISO2 => begin move_to_iiso(x); infiso2 += 1; end
                :HOS2  => begin move_to_hospicu(x); hos2 += 1; end 
                :ICU2  => begin move_to_hospicu(x); icu2 += 1; end
                :REC2  => begin move_to_recovered(x); rec2 += 1; end
                :DED2  => begin move_to_dead(x); ded2 += 1; end
                :LAT3  => begin move_to_latent(x); lat3 += 1; end
                :PRE3  => begin move_to_pre(x); pre3 += 1; end
                :ASYMP3 => begin move_to_asymp(x); asymp3 += 1; end
                :MILD3 => begin move_to_mild(x); mild3 += 1; end
                :MISO3 => begin move_to_miso(x); miso3 += 1; end
                :INF3  => begin move_to_inf(x); inf3 +=1; end    
                :IISO3 => begin move_to_iiso(x); infiso3 += 1; end
                :HOS3  => begin move_to_hospicu(x); hos3 += 1; end 
                :ICU3  => begin move_to_hospicu(x); icu3 += 1; end
                :REC3  => begin move_to_recovered(x); rec3 += 1; end
                :DED3  => begin move_to_dead(x); ded3 += 1; end
                _    => error("swap expired, but no swap set.")
            end
        end
        # run covid-19 functions for other integrated dynamics. 
        #ct_dynamics(x)

        # get the meet counts for the next day 
        get_nextday_counts(x)
        if p.vaccinating
            vac_update(x)
        end
    end
    return (lat, mild, miso, inf, infiso, hos, icu, rec, ded,lat2, mild2, miso2, inf2, infiso2, hos2, icu2, rec2, ded2,lat3, mild3, miso3, inf3, infiso3, hos3, icu3, rec3, ded3)
end
export time_update

@inline _set_isolation(x::Human, iso) = _set_isolation(x, iso, x.isovia)
@inline function _set_isolation(x::Human, iso, via)
    # a helper setter function to not overwrite the isovia property. 
    # a person could be isolated in susceptible/latent phase through contact tracing
    # --> in which case it will follow through the natural history of disease 
    # --> if the person remains susceptible, then iso = off
    # a person could be isolated in presymptomatic phase through fpreiso
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # a person could be isolated in mild/severe phase through fmild, fsevere
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # --> if x.iso == true from PRE and x.isovia == :pi, do not overwrite
    x.iso = iso 
    x.isovia == :null && (x.isovia = via)
end

function sample_epi_durations()
    # when a person is sick, samples the 
    #lat_dist = Distributions.truncated(Gamma(3.122, 2.656),4,11.04) # truncated between 4 and 7
    lat_dist = Distributions.truncated(LogNormal(1.434, 0.661), 4, 7) # truncated between 4 and 7
    pre_dist = Distributions.truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = Int.(round.(rand(lat_dist)))
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = Int.(ceil.(rand(asy_dist)))
    infs = Int.(ceil.(rand(inf_dist)))
    return (latents, asymps, pres, infs)
end

function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration
    x.health = x.swap
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
   
    #0-18 31 19 - 59 29 60+ 18 going to asymp
    symp_pcts = [0.7, 0.623, 0.672, 0.672, 0.812, 0.812] #[0.3 0.377 0.328 0.328 0.188 0.188]
    age_thres = [4, 19, 49, 64, 79, 999]
    g = findfirst(y-> y >= x.age, age_thres)

    if x.recovered
        auxiliar = (1-p.vac_efficacy_symp[1][2][end])
    else
        aux = x.vac_status*x.protected > 0 ? p.vac_efficacy_symp[x.strain][x.vac_status][x.protected] : 0.0
        auxiliar = (1-aux)
    end
 
    if rand() < (symp_pcts[g])*auxiliar
        if x.strain == 1
            x.swap = PRE
        elseif x.strain == 2
            x.swap = PRE2
        elseif x.strain == 3
            x.swap = PRE3
        else
            error("No strain in move to lat")
        end
        #x.swap = x.strain == 1 ? PRE : PRE2
    else
        if x.strain == 1
            x.swap = ASYMP
        elseif x.strain == 2
            x.swap = ASYMP2
        elseif x.strain == 3
            x.swap = ASYMP3
        else
            error("No strain in move to lat")
        end

        #x.swap = x.strain == 1 ? ASYMP : ASYMP2
    end
    x.wentTo = x.swap
    x.got_inf = true
    ## in calibration mode, latent people never become infectious.
    if p.calibration && !x.first_one
        x.swap = LAT 
        x.exp = 999
    end 
end
export move_to_latent

function move_to_asymp(x::Human)
    ## transfers human h to the asymptomatic stage 
    x.health = x.swap  
    x.tis = 0 
    x.exp = x.dur[2] # get the presymptomatic period
   

    if x.strain == 1
        x.swap = REC
    elseif x.strain == 2
        x.swap = REC2
    elseif x.strain == 3
        x.swap = REC3
    else
        error("No strain in move to asymp")
    end
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, the asymptomatic individual has limited contacts
end
export move_to_asymp

function move_to_pre(x::Human)
    θ = (0.95, 0.9, 0.85, 0.6, 0.2)  # percentage of sick individuals going to mild infection stage
    x.health = x.swap
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period

    if x.recovered
        auxiliar = (1-p.vac_efficacy_sev[1][2][end])
    else
        aux = x.vac_status*x.protected > 0 ? p.vac_efficacy_sev[x.strain][x.vac_status][x.protected] : 0.0
        auxiliar = (1-aux)
    end

    if rand() < (1-θ[x.ag])*auxiliar
        if x.strain == 1
            x.swap = INF
        elseif x.strain == 2
            x.swap = INF2
        elseif x.strain == 3
            x.swap = INF3
        else
            error("No strain in move to pre")
        end
        #x.swap = x.strain == 1 ? INF : INF2
    else 
        if x.strain == 1
            x.swap = MILD
        elseif x.strain == 2
            x.swap = MILD2
        elseif x.strain == 3
            x.swap = MILD3
        else
            error("No strain in move to pre")
        end
        #x.swap = x.strain == 1 ? MILD : MILD2
    end
    # calculate whether person is isolated
    rand() < p.fpreiso && _set_isolation(x, true, :pi)
end
export move_to_pre

function move_to_mild(x::Human)
    ## transfers human h to the mild infection stage for γ days
   
    x.health = x.swap 
    x.tis = 0 
    x.exp = x.dur[4]
    if x.strain == 1
        x.swap = REC
    elseif x.strain == 2
        x.swap = REC2
    elseif x.strain == 3
        x.swap = REC3
    else
        error("No strain in move to mild")
    end
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, staying in MILD is same as MISO since contacts will be limited. 
    # we still need the separation of MILD, MISO because if x.iso is false, then here we have to determine 
    # how many days as full contacts before self-isolation
    # NOTE: if need to count non-isolated mild people, this is overestimate as isolated people should really be in MISO all the time
    #   and not go through the mild compartment 
    aux = x.vac_status > 0 ? p.fmild*p.red_risk_perc : p.fmild
    if x.iso || rand() < aux#p.fmild
        if x.strain == 1
            x.swap = MISO
        elseif x.strain == 2
            x.swap = MISO2
        elseif x.strain == 3
            x.swap = MISO3
        else
            error("No strain in move to mild")
        end
        #x.swap = x.strain == 1 ? MISO : MISO2  
        x.exp = p.τmild
    end
end
export move_to_mild

function move_to_miso(x::Human)
    ## transfers human h to the mild isolated infection stage for γ days
    x.health = x.swap
    if x.strain == 1
        x.swap = REC 
    elseif x.strain == 2
        x.swap = REC2
    elseif x.strain == 3
        x.swap = REC3 
    else
        error("No strain in move to miso")
    end
    #x.swap = x.strain == 1 ? REC : REC2
    x.tis = 0 
    x.exp = x.dur[4] - p.τmild  ## since tau amount of days was already spent as infectious
    _set_isolation(x, true, :mi) 
end
export move_to_miso

function move_to_infsimple(x::Human)
    ## transfers human h to the severe infection stage for γ days 
    ## simplified function for calibration/general purposes
    x.health = x.swap
    x.tis = 0 
    x.exp = x.dur[4]
    if x.strain == 1
        x.swap = REC 
    elseif x.strain == 2
        x.swap = REC2
    elseif x.strain == 3
        x.swap = REC3 
    else
        error("No strain in move to miso")
    end
    #x.swap = x.strain == 1 ? REC : REC2
    _set_isolation(x, false, :null) 
end

function move_to_inf(x::Human)
    ## transfers human h to the severe infection stage for γ days
    ## for swap, check if person will be hospitalized, selfiso, die, or recover
 
    # h = prob of hospital, c = prob of icu AFTER hospital    
   
    h = x.comorbidity == 1 ? 1.0 : 0.09 #0.376
    c = x.comorbidity == 1 ? 0.33 : 0.25
    
    groups = [0:34,35:54,55:69,70:84,85:100]
    gg = findfirst(y-> x.age in y,groups)

    mh = [0.0002; 0.0015; 0.011; 0.0802; 0.381] # death rate for severe cases.
   
    ###prop/(prob de sintoma severo)
    if p.calibration && !p.calibration2
        h =  0#, 0, 0, 0)
        c =  0#, 0, 0, 0)
        mh = (0, 0, 0, 0, 0)
    end

    time_to_hospital = Int(round(rand(Uniform(2, 5)))) # duration symptom onset to hospitalization
   	
    x.health = x.swap
    x.swap = UNDEF
    x.tis = 0 
    if rand() < h     # going to hospital or ICU but will spend delta time transmissing the disease with full contacts 
        x.exp = time_to_hospital
        if rand() < c
            if x.strain == 1
                x.swap = ICU
            elseif x.strain == 2
                x.swap = ICU2
            elseif x.strain == 3
                x.swap = ICU3
            else
                error("No strain in move to inf")
            end
            #x.swap = x.strain == 1 ? ICU : ICU2
        else
            if x.strain == 1
                x.swap = HOS
            elseif x.strain == 2
                x.swap = HOS2
            elseif x.strain == 3
                x.swap = HOS3
            else
                error("No strain in move to inf")
            end
            #x.swap = x.strain == 1 ? HOS : HOS2
        end
       
    else ## no hospital for this lucky (but severe) individual 
        aux = (p.mortality_inc^Int(x.strain==2))
        
        if x.iso || rand() < p.fsevere 
            x.exp = 1  ## 1 day isolation for severe cases 
            if x.strain == 1
                x.swap = IISO
            elseif x.strain == 2
                x.swap = IISO2
            elseif x.strain == 3
                x.swap = IISO3
            else
                error("No strain in move to inf")
            end    
            #x.swap = x.strain == 1 ? IISO : IISO2
        else
            if rand() < mh[gg]*aux
                x.exp = x.dur[4] 
                if x.strain == 1
                    x.swap = DED
                elseif x.strain == 2
                    x.swap = DED2
                elseif x.strain == 3
                    x.swap = DED3
                else
                    error("No strain in move to inf")
                end 
            else 
                x.exp = x.dur[4]  
                if x.strain == 1
                    x.swap = REC 
                elseif x.strain == 2
                    x.swap = REC2
                elseif x.strain == 3
                    x.swap = REC3 
                else
                    error("No strain in move to miso")
                end
            end

        end  
       
    end
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent I -> ?")
end

function move_to_iiso(x::Human)
    ## transfers human h to the sever isolated infection stage for γ days
    x.health = x.swap
    groups = [0:34,35:54,55:69,70:84,85:100]
    gg = findfirst(y-> x.age in y,groups)
    mh = [0.0002; 0.0015; 0.011; 0.0802; 0.381] # death rate for severe cases.
    aux = (p.mortality_inc^Int(x.strain==2))
    if rand() < mh[gg]*aux
        x.exp = x.dur[4] 
        if x.strain == 1
            x.swap = DED
        elseif x.strain == 2
            x.swap = DED2
        elseif x.strain == 3
            x.swap = DED3
        else
            error("No strain in move to inf")
        end 
    else 
        x.exp = x.dur[4]  
        if x.strain == 1
            x.swap = REC 
        elseif x.strain == 2
            x.swap = REC2
        elseif x.strain == 3
            x.swap = REC3 
        else
            error("No strain in move to miso")
        end
    end
    #x.swap = x.strain == 1 ? REC : REC2
    x.tis = 0     ## reset time in state 
    x.exp = x.dur[4] - 1  ## since 1 day was spent as infectious
    _set_isolation(x, true, :mi)
end 

function move_to_hospicu(x::Human)   
    #death prob taken from https://www.cdc.gov/nchs/nvss/vsrr/covid_weekly/index.htm#Comorbidities
    # on May 31th, 2020
    #= age_thres = [24;34;44;54;64;74;84;999]
    g = findfirst(y-> y >= x.age,age_thres) =#
    aux = [0:4, 5:19, 20:44, 45:54, 55:64, 65:74, 75:84, 85:99]
   
    if x.strain == 1 || x.strain == 3

        mh = [0.001, 0.001, 0.0015, 0.0065, 0.01, 0.02, 0.0735, 0.38]
        mc = [0.002,0.002,0.0022, 0.008, 0.022, 0.04, 0.08, 0.4]

    elseif x.strain == 2
    
        mh = [0.001, 0.001, 0.0025, 0.008, 0.02, 0.038, 0.15, 0.66]
        mc = [0.002,0.002,0.0032, 0.01, 0.022, 0.04, 0.2, 0.70]

    else
      
            error("No strain - hospicu")
    end
    
    gg = findfirst(y-> x.age in y,aux)

    psiH = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17))))
    psiC = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17)))) + 2
    muH = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15))))
    muC = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15)))) + 2

    swaphealth = x.swap 
    x.health = swaphealth ## swap either to HOS or ICU
    x.swap = UNDEF
    x.tis = 0
    _set_isolation(x, true) # do not set the isovia property here.  

    if swaphealth == HOS || swaphealth == HOS2 || swaphealth == HOS3
        x.hospicu = 1 
        if rand() < mh[gg] ## person will die in the hospital 
            x.exp = muH 
            if x.strain == 1
                x.swap = DED
            elseif x.strain == 2
                x.swap = DED2
            elseif x.strain == 3
                x.swap = DED3
            else
                error("No strain in move to inf")
            end 
            #x.swap = x.strain == 1 ? DED : DED2
        else 
            x.exp = psiH 
            if x.strain == 1
                x.swap = REC 
            elseif x.strain == 2
                x.swap = REC2
            elseif x.strain == 3
                x.swap = REC3 
            else
                error("No strain in move to miso")
            end
            #x.swap = x.strain == 1 ? REC : REC2
        end    
    elseif swaphealth == ICU || swaphealth == ICU2 || swaphealth == ICU3
        x.hospicu = 2 
                
        if rand() < mc[gg] ## person will die in the ICU 
            x.exp = muC
            if x.strain == 1
                x.swap = DED
            elseif x.strain == 2
                x.swap = DED2
            elseif x.strain == 3
                x.swap = DED3
            else
                error("No strain in move to inf")
            end 
            #x.swap = x.strain == 1 ? DED : DED2
        else 
            x.exp = psiC
            if x.strain == 1
                x.swap = REC 
            elseif x.strain == 2
                x.swap = REC2
            elseif x.strain == 3
                x.swap = REC3 
            else
                error("No strain in move to miso")
            end
            #x.swap = x.strain == 1 ? REC : REC2
        end
    else
        error("error in hosp")
    end
    
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent H -> ?")    
end

function move_to_dead(h::Human)
    # no level of alchemy will bring someone back to life. 
    h.health = h.swap
    h.swap = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.iso = true # a dead person is isolated
    _set_isolation(h, true)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end

function move_to_recovered(h::Human)
    h.health = h.swap

    if h.strain in (1,2)
        h.recovered = true
    end

    h.swap = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.iso = false ## a recovered person has ability to meet others
    _set_isolation(h, false)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end


@inline function _get_betavalue(sys_time, xhealth) 
    #bf = p.β ## baseline PRE
    length(BETAS) == 0 && return 0
    bf = BETAS[sys_time]
    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
    if xhealth == ASYMP
        bf = bf * p.frelasymp #0.11

    elseif xhealth == MILD || xhealth == MISO 
        bf = bf * 0.44

    elseif xhealth == INF || xhealth == IISO 
        bf = bf * 0.89

    elseif xhealth == ASYMP2
        bf = bf*p.frelasymp*p.sec_strain_trans #0.11

    elseif xhealth == MILD2 || xhealth == MISO2
        bf = bf * 0.44*p.sec_strain_trans

    elseif xhealth == INF2 || xhealth == IISO2 
        bf = bf * 0.89*p.sec_strain_trans

    elseif xhealth == PRE2
        bf = bf*p.sec_strain_trans

    elseif xhealth == ASYMP3
        bf = bf*p.frelasymp*p.third_strain_trans #0.11

    elseif xhealth == MILD3 || xhealth == MISO3
        bf = bf * 0.44*p.third_strain_trans

    elseif xhealth == INF3 || xhealth == IISO3 
        bf = bf * 0.89*p.third_strain_trans

    elseif xhealth == PRE3
        bf = bf*p.third_strain_trans
    end
    return bf
end
export _get_betavalue

@inline function get_nextday_counts(x::Human)
    # get all people to meet and their daily contacts to recieve
    # we can sample this at the start of the simulation to avoid everyday    
    cnt = 0
    ag = x.ag
    #if person is isolated, they can recieve only 3 maximum contacts
    
    if !x.iso 
        #cnt = rand() < 0.5 ? 0 : rand(1:3)
        aux = x.relaxed ? 1.0 : p.contact_change_rate*p.contact_change_2
        cnt = rand(negative_binomials(ag,aux)) ##using the contact average for shelter-in
    else 
        cnt = rand(negative_binomials_shelter(ag,p.contact_change_2))  # expensive operation, try to optimize
    end
    
    if x.health in (DED,DED2,DED3)
        cnt = 0 
    end
    x.nextday_meetcnt = cnt
    return cnt
end

function dyntrans(sys_time, grps,sim)
    totalmet = 0 # count the total number of contacts (total for day, for all INF contacts)
    totalinf = 0 # count number of new infected 
    ## find all the people infectious
    #rng = MersenneTwister(246*sys_time*sim)
    pos = shuffle(1:length(humans))
    # go through every infectious person
    for x in humans[pos]        
        if x.health in (PRE, ASYMP, MILD, MISO, INF, IISO,PRE2, ASYMP2, MILD2, MISO2, INF2, IISO2, PRE3, ASYMP3, MILD3, MISO3, INF3, IISO3)
            
            xhealth = x.health
            cnts = x.nextday_meetcnt
            cnts == 0 && continue # skip person if no contacts
            
            gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
            for (i, g) in enumerate(gpw) 
                meet = rand(grps[i], g)   # sample the people from each group
                # go through each person
                for j in meet 
                    y = humans[j]
                    ycnt = y.nextday_meetcnt    
                    ycnt == 0 && continue

                    y.nextday_meetcnt = y.nextday_meetcnt - 1 # remove a contact
                    totalmet += 1
                    
                    beta = _get_betavalue(sys_time, xhealth)
                    adj_beta = 0 # adjusted beta value by strain and vaccine efficacy
                    if y.health == SUS && y.swap == UNDEF                  
                        aux = y.vac_status*y.protected > 0 ? p.vac_efficacy_inf[x.strain][y.vac_status][y.protected] : 0.0
                        adj_beta = beta*(1-aux)
                    elseif (x.strain == 3 && y.health in (REC, REC2) && y.swap == UNDEF)
                        adj_beta = beta*(p.reduction_recovered) #0.21
                    end

                    if rand() < adj_beta
                        totalinf += 1
                        y.exp = y.tis   ## force the move to latent in the next time step.
                        y.sickfrom = xhealth ## stores the infector's status to the infectee's sickfrom
                        y.sickby = x.idx
                        y.strain = x.strain       
                        if x.strain == 1
                            y.swap = LAT
                        elseif x.strain == 2
                            y.swap = LAT2
                        elseif x.strain == 3 
                            y.swap = LAT3 
                        else
                            error("No strain in move to transmission")
                        end 
                        #y.swap = y.strain == 1 ? LAT : LAT2
                    end  
                end
            end            
        end
    end
    return totalmet, totalinf
end
export dyntrans

function contact_matrix()
    # regular contacts, just with 5 age groups. 
    #  0-4, 5-19, 20-49, 50-64, 65+
    CM = Array{Array{Float64, 1}, 1}(undef, 5)
     CM[1] = [0.2287, 0.1839, 0.4219, 0.1116, 0.0539]
    CM[2] = [0.0276, 0.5964, 0.2878, 0.0591, 0.0291]
    CM[3] = [0.0376, 0.1454, 0.6253, 0.1423, 0.0494]
    CM[4] = [0.0242, 0.1094, 0.4867, 0.2723, 0.1074]
    CM[5] = [0.0207, 0.1083, 0.4071, 0.2193, 0.2446] 
   
    return CM
end
# 
# calibrate for 2.7 r0
# 20% selfisolation, tau 1 and 2.

function negative_binomials(ag,mult) 
    ## the means/sd here are calculated using _calc_avgag
    means = [10.21, 16.793, 13.7950, 11.2669, 8.0027]
    sd = [7.65, 11.7201, 10.5045, 9.5935, 6.9638]
    means = means*mult
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]
end
#const nbs = negative_binomials()
const cm = contact_matrix()
#export negative_binomials, contact_matrix, nbs, cm

export negative_binomials


function negative_binomials_shelter(ag,mult) 
    ## the means/sd here are calculated using _calc_avgag
    means = [2.86, 4.7, 3.86, 3.15, 2.24]
    sd = [2.14, 3.28, 2.94, 2.66, 1.95]
    means = means*mult
    #sd = sd*mult
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]   
end

## internal functions to do intermediate calculations
function _calc_avgag(lb, hb) 
    ## internal function to calculate the mean/sd of the negative binomials
    ## returns a vector of sampled number of contacts between age group lb to age group hb
    dists = _negative_binomials_15ag()[lb:hb]
    totalcon = Vector{Int64}(undef, 0)
    for d in dists 
        append!(totalcon, rand(d, 10000))
    end    
    return totalcon
end
export _calc_avgag

function _negative_binomials_15ag()
    ## negative binomials 15 agegroups
    AgeMean = Vector{Float64}(undef, 15)
    AgeSD = Vector{Float64}(undef, 15)
    #0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70+
    #= AgeMean = [10.21, 14.81, 18.22, 17.58, 13.57, 13.57, 14.14, 14.14, 13.83, 13.83, 12.3, 12.3, 9.21, 9.21, 6.89]
    AgeSD = [7.65, 10.09, 12.27, 12.03, 10.6, 10.6, 10.15, 10.15, 10.86, 10.86, 10.23, 10.23, 7.96, 7.96, 5.83]
     =#
     AgeMean = repeat([14.14],15)#[10.21, 14.81, 18.22, 17.58, 13.57, 13.57, 14.14, 14.14, 13.83, 13.83, 12.3, 12.3, 9.21, 9.21, 6.89]
    AgeSD = repeat([10.86],15)#[7.65, 10.09, 12.27, 12.03, 10.6, 10.6, 10.15, 10.15, 10.86, 10.86, 10.23, 10.23, 7.96, 7.96, 5.83]
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, 15)
    for i = 1:15
        p = 1 - (AgeSD[i]^2-AgeMean[i])/(AgeSD[i]^2)
        r = AgeMean[i]^2/(AgeSD[i]^2-AgeMean[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms    
end



const vaccination_days = days_vac_f()
const vac_rate_1 = vaccination_rate_1()
const vac_rate_2 = vaccination_rate_2()
## references: 
# critical care capacity in Canada https://www.ncbi.nlm.nih.gov/pubmed/25888116
end # module end
