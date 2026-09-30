#packages necessaires
library(nnet)
library(tidyverse)
library(stringr)


###########################################################################################################
#####################                 BuildMMEstInputs                #####################################
###########################################################################################################

CheckNamesModele<-function(data,fixed,random,error){
  
  #Check if fixed is a character
  if(!is.character(fixed)){
    stop("Fixed should be a character")
  }
  
  #Check if random is a character
  if(!is.character(random)){
    stop("Random should be a character")
  }
  
  #Check if error is a character
  if(!is.character(error)){
    stop("Error should be a character")
  }
  
  
  # No colnames name "error" 
  if("error" %in% colnames(data)){
    stop("You can't use 'error' as colnames in data")
  }
  
  # No colnames name "Xeffect" 
  if("Xeffect" %in% colnames(data)){
    stop("You can't use 'Xeffect' as colnames in data")
  }
  
  # compilation of all formulas
  total<-paste(fixed,random,error)
  
  #names cannot include '+'. '+' should be between two white space.
  if(!all(sapply(str_locate_all(total, "[[+]]")[[1]][,1],function(x){str_sub(total,x-1,x+1)} == " + "))){
    stop("Please, Use a space between effect and '+' in fixed, random and error or remove '+' in your names")
  }
  
  #names cannot include '~'. '~' should be between two white space.
  if(!all(sapply(str_locate_all(fixed, "~")[[1]][,1],function(x){str_sub(total,x-1,x+1)} == " ~ "))){
    stop("Please, Use a space between your response and "~" in fixed or remove '~' in your names")
  }
  
  #formula should start with a '~' except for fixed.
  if(!all(unlist(lapply(list(random,error),function(x){substring(x, 1, 1) == "~"})))){
    stop("Please, Begin your fomula in 'random' and 'error' by '~'")
  }
  
  #'~' should not be used in formulas, except at the beginning.
  if(any(unlist(lapply(list(random,error),function(x){str_detect(x,"~") & substring(x, 1, 1) != "~"})))){
    stop("Unappropriate use of '~' in random or error")
  }
  
  # random or error do not handle interaction with X
  if(any(unlist(lapply(list(random,error),function(x){str_detect(x,"Xeffect")})))){
    stop("Random or error do not handle interaction with X")
  }
  
  #Which character chain are forbidden ?
  forbidden<-c("at(","vm(","~mu",":","(",")","_x_","Xeffect")
  Names<-str_remove_all(total,c(" ")) %>% str_remove_all(.,c("error"))  %>%
    str_split(.,"~mu|~at\\(vm\\(|~vm\\(at\\(|~at\\(|~vm\\(|~|[[+]]at\\(vm\\(|[[+]]vm\\(at\\(|[[+]]at\\(|[[+]]vm\\(|[[+]]|\\)|:|Xeffect") %>% .[[1]] %>% .[-which(. == "")]
  check_forbidden<-sapply(forbidden,function(x){str_detect(Names,fixed(x))}) %>% matrix(.,ncol=length(forbidden)) %>%  apply(.,2,any)
  if(any(check_forbidden)){
    stop(paste0("Unappropriate use of ", names(check_forbidden[check_forbidden == T]), "\n"))  
  }
  
}

#Function for checking if effect names exist in data
CheckColnames<-function(data,ColNames){
  if(!any(ColNames%in%colnames(data))){
    stop(paste0("Wrong name for ",ColNames))
  }
}

#function to calculate each indicator for fixed effect without interactions
CalcSimpleMatrixFixed<-function(data,Effect){
  
  #saving of effect names
  name_effect<-Effect
  
  #Extraction of the data vector of the fixed effect
  Effect<-data[,Effect]
  
  #check on the class of the effect : should be numeric or factor
  if(class(Effect)=="character"){
    Effect<-as.factor(Effect)
  }
  if(class(Effect)%in% c("numeric","integer","double")){
    Effect<-as.numeric(Effect)
  }
  
  #If it's a factor, we build on indicator for aech modality
  if(class(Effect)=="factor"){
    Effect<-factor(Effect,levels=unique(Effect))
    Effect<-class.ind(Effect)
  }else{#if it is numeric, we use data vector as an indicator
    Effect<-Effect
  }
  
  #Formating
  Effect<-as.data.frame(Effect)
  if(ncol(Effect)==1){
    colnames(Effect)<-name_effect
  }
  
  return(Effect)
}

#Function to convert complicated names for interaction effect into a list easy to use
FormatingInteraction<-function(x){
  
  #x<-interaction[1]
  
  x<-unlist(str_split(x,":"))%>%
    str_split(.,",")
  
  x[[1]]<-str_remove_all(x[[1]],"at\\(|\\)")
  if(length(x[[1]])==1){#all level
    names(x)<-c("all_level","principal_effect")
    names(x[[1]])<-"secondary_effect"
  }else{#specific level
    names(x)<-c("specific_level","principal_effect")
    names(x[[1]])<-c("secondary_effect","level")
  }
  return(x)
  
}

#function to construct indicators for fixed effect with interaction
CalcInteractionMatrixFixed<-function(Effect_interaction,Liste_fixed){
  
  #Effect_interaction<-interaction[[3]]
  
  #To avoid useless calculation, we check if the effect indicator is already build. 
  if(Effect_interaction$principal_effect %in% names(Liste_fixed)){
    principal <- Liste_fixed[[Effect_interaction$principal_effect]]%>%
      as.matrix(.,drop=FALSE)
    
  }else{#else, we need to use the CalcSimpleMatrixFixed to calculated the indicator effect
    principal<-lapply(Effect_interaction$principal_effect,function(x){CalcSimpleMatrixFixed(data,x)}) %>%
      map2(.x=.,.y=Effect_interaction$principal_effect, ~ {colnames(.x)<-paste0(.y,"_x_",colnames(.x));return(.x)})%>%
      .[[1]]%>%
      as.matrix(.,drop=FALSE)
  }
  
  # Interaction effect can be global or conditional to specific modality
  
  if(names(Effect_interaction)[1]=="all_level"){ # interaction with all levels (Global)
    
    if(Effect_interaction$all_level[1] %in% names(Liste_fixed)){#To avoid useless calculation, we check if the effect indicator is already build. 
      secondary <- Liste_fixed[[Effect_interaction$all_level[1]]]%>%
        as.matrix(.,drop=FALSE)
    }else{#else, we need to use the CalcSimpleMatrixFixed to calculated the indicator effect
      secondary<-lapply(Effect_interaction$all_level[1],function(x){CalcSimpleMatrixFixed(data,x)}) %>%
        map2(.x=.,.y=Effect_interaction$all_level[1], ~ {colnames(.x)<-paste0(.y,"_x_",colnames(.x));return(.x)})%>%
        .[[1]]%>%
        as.matrix(.,drop=FALSE)
    }
  }
  
  
  if(names(Effect_interaction)[1]=="specific_level"){# interacion with a specific level
    
    if(Effect_interaction$specific_level[1] %in% names(Liste_fixed)){#To avoid useless calculation, we check if the effect indicator is already build.
      secondary <- Liste_fixed[[Effect_interaction$specific_level[1]]]%>%
        as.matrix(.,drop=FALSE)
    }else{#else, we need to use the CalcSimpleMatrixFixed to calculated the indicator effect
      secondary<-lapply(Effect_interaction$specific_level[1],function(x){CalcSimpleMatrixFixed(data,x)}) %>%
        map2(.x=.,.y=Effect_interaction$specific_level[1], ~ {colnames(.x)<-paste0(.y,"_x_",colnames(.x));return(.x)})%>%
        .[[1]]%>%
        as.matrix(.,drop=FALSE)
    }
    secondary<-secondary[,paste0(Effect_interaction$specific_level[1],"_x_",Effect_interaction$specific_level[2]),drop=FALSE]
  }
  
  #construction of interaction indicators
  matrix_interaction<-lapply(1:ncol(principal),function(y){apply(secondary,2,function(x){principal[,y] * x})})%>%
    do.call("cbind",.)
  colnames(matrix_interaction)<-paste0(rep(colnames(principal),each=ncol(secondary)),":",colnames(secondary))
  
  return(matrix_interaction)
}

#function to construct design matrix and covariance matrix for random effect without interaction
CalcSimpleMatrixRandom<-function(data,x){
  
  #check on effect name
  CheckColnames(data,x[1])
  
  #check on the class of the effect : should be factor
  if(class(data[,x[1]])=="character"){
    data[,x[1]]<-as.factor(data[,x[1]])
  }
  if(class(data[,x[1]])!="factor"){
    stop(paste0(x[1])," should be a factor or a character")
  }
  
  #Formating factor
  data[,x[1]]<-factor(data[,x[1]],levels=unique(data[,x[1]]))
  
  #contruction of the design matrix
  Z<-class.ind(data[,x[1]])
  
  if(length(x)==1){#If no covarance matrix is defined
    V<-diag(ncol(Z)) #the covariance matrix is the diagonal matrix
    dimnames(V)<-list(colnames(Z),colnames(Z))
  }else{#If the covariance matrix is defined
    if(exists(x[2])){ #If the covariance matrix exist in the global environnement
      K<-get(x[2],envir = .GlobalEnv)
    }else{
      stop(paste0(x[2]," does not exist"))
    }
    if(!all(colnames(Z) %in% colnames(K))){ #check of the adequation of colnames between K and Z
      #stop(paste0(x[2]," should have the same names/levels as ",x[1]))
      warning(paste0(x[2]," has not the same names/levels as ",x[1]))
      Z<-Z[,which(colnames(Z) %in% colnames(K))]
    }
    V<-K[colnames(Z),colnames(Z)] #order K following Z columns
  }
  
  List_matrix<-list(Z,V)
  names(List_matrix)<-c("Z","V")
  return(List_matrix)
}

#function to construct design matrix and covariance matrix for random effect with interaction
CalcInteractionMatrixRandom<-function(data,x,Liste_random,Liste_fixed){
  
  if(x[2] %in% str_remove_all(names(Liste_random),pattern="\\(|\\)") ){ #check if the effect is already done
    Z_random<-Liste_random[[paste0("(",x[2],")")]][["Z"]]
    
  }else{#Else we construct the design matrice and the covariance matrix using CalcSimpleMatrixRandom function
    temp<-CalcSimpleMatrixRandom(data,str_split(x[2],",")[[1]])
    Z_random<-temp[[1]]
  }
  
  K<-str_split(x[2],",")[[1]][2]
  if(is.na(K)){ #if K is not defined
    K<-diag(ncol(Z_random)) #K is the diagonal matrix
    dimnames(K)<-list(colnames(Z_random),colnames(Z_random))
  }else{#else, we extract K from the global environnement
    if(exists(K)){ #If the covariance matrix exist in the global environnement
      K<-get(K,envir = .GlobalEnv)
    }else{
      stop(paste0(K,"doesnt exist"))
    }
  }
  
  if(x[1] %in% names(Liste_fixed)){#If we have already build the fixed effect part of the interaction
    Z_fixed<-Liste_fixed[[x[1]]]
  }else{#else we need to build it with CalcSimpleMatrixFixed function
    temp<-str_remove_all(x[1],"at\\(") %>% str_split(.,",") %>%.[[1]]
    Z_fixed<-CalcSimpleMatrixFixed(data,temp[1])
    colnames(Z_fixed)<-paste0(temp[1],"_x_",colnames(Z_fixed))
    if(length(temp)!=1){
      Z_fixed<-Z_fixed[,paste0(temp[1],"_x_",temp[2]),drop=FALSE]
    }
  }
  
  #construction of the interaction design matrix
  Z_specific<-lapply(lapply(1:ncol(Z_fixed),function(x){rep(x,dim(Z_random)[2])}),function(x){Z_fixed[,x]}) %>%
    map(.,~.*Z_random) %>%
    lapply(.,function(x){colnames(x)<-colnames(Z_random);x})
  names(Z_specific)<-colnames(Z_fixed)
  
  if(!str_detect(x[1],"at\\(")){# the same variance to each modality
    
    Z<-do.call("cbind",Z_specific)
    colnames(Z)<-str_replace_all(colnames(Z),"[[.]]",paste0(":",str_split(x[2],",")[[1]][1],"_x_"))
    V<-kronecker(diag(length(Z_specific)),K[colnames(Z_random),colnames(Z_random)])
    dimnames(V)<-list(colnames(Z),colnames(Z))
    
    Z<-list(Z)
    names(Z)<-str_remove_all(x[1],"at\\(") %>% str_split(.,",") %>%.[[1]]
    V<-list(V)
    names(V)<-str_remove_all(x[1],"at\\(") %>% str_split(.,",") %>%.[[1]]
    
    
  }else{# specific variance to each modality
    
    Z<-Z_specific
    V<-lapply(1:length(Z_specific),function(x){return(K[colnames(Z_random),colnames(Z_random)])})
    Z<-imap(Z,function(.x,.y){colnames(.x)<-paste0(.y,":",str_split(x[2],",")[[1]][1],"_x_",colnames(.x));return(.x)})
    V<-map2(Z,V,function(.x,.y){dimnames(.y)<-list(colnames(.x),colnames(.x));return(.y)})
  }
  
  List_matrix<-list(lapply(Z,as.matrix),lapply(V,as.matrix))
  names(List_matrix)<-c("Z","V")
  
  return(List_matrix)
}

#function to compilate all random effet matrix. Can work for the design matrix or the covariance matrix.
compilation_random<-function(Liste_random,Liste_random_interaction,interaction,Type){
  
  # Random effects with no interaction
  if(!is.null(Liste_random)){
    liste_random<-lapply(1:length(Liste_random),function(x){Liste_random[[x]][[Type]]})
    names(liste_random)<-names(Liste_random)
  }else{
    liste_random<-NULL
  }
  
  # Random effects with interaction
  if(!is.null(Liste_random_interaction)){
    list_random_interaction<-lapply(1:length(Liste_random_interaction),function(x){Liste_random_interaction[[x]][[Type]]})
    names(list_random_interaction)<-names(Liste_random_interaction)<-lapply(interaction,function(x){x[2]}) %>%
      paste0("(",.,")")
    list_random_interaction<-unlist(list_random_interaction, recursive=FALSE)
    names(list_random_interaction)<-str_replace_all(names(list_random_interaction),"[[.]]",":")
  }else{
    list_random_interaction<-NULL
  }
  
  #compilation
  if(!is.null(liste_random)){
    if(!is.null(list_random_interaction)){
      List<-c(liste_random,list_random_interaction)
    }else{
      List<-c(liste_random)
    }
  }else{
    if(!is.null(list_random_interaction)){
      List<-c(list_random_interaction)
    }else{
      List<-NULL
    }
  }
  
  return(List)
}

BuildMMEstInputs<-function(data,X=NULL,fixed,random="~",error="~",CheckNames=FALSE){
  
  #####################################################
  ##############          Check          ##############
  #####################################################
  
  ##### Check on data names use in the modele
  
  if(CheckNames){
    CheckNamesModele(data,fixed,random,error)
  }
  
  if(!is.null(X)){
    if(!is.matrix(X)){
      stop("X should be a matrix")
    }
  }
  
  
  ##### Check on fixed effects (use in NA treatment)
  
  #Xeffect should not be in fixed with any X matrix
  if(is.null(X) & str_detect(fixed,"Xeffect")){
    stop("Xeffect is specified in 'fixed' with any X matrix provided")
  }
  
  #extraction of fixed effects
  fixed_effect<-str_split(fixed,"[[~]]")[[1]][2] %>%
    str_replace_all(., " ", "") %>%
    str_split(.,"[[+]]")%>%
    .[[1]]
  
  fixed_effect<-fixed_effect[which(!str_detect(fixed_effect,"Xeffect"))]
  
  #check of names effect
  if(all(fixed_effect != "")){
    fixed_effect_split<-str_split(fixed_effect,":") %>% unlist %>% unique
    lapply(fixed_effect_split,function(x){CheckColnames(data,x)})
  }
  
  #####################################################
  ##############      NA treatment       ##############
  #####################################################
  
  
  line_noNA_Y<-which(!is.na(data[,str_split(fixed,"[[~]]")[[1]][1] %>% str_replace_all(., " ", "")]))
  NA_Y<-nrow(data)-length(line_noNA_Y)
  
  if(NA_Y != 0){
    warning(NA_Y," observations containing missing values in the response are elimnated.")
  }
  
  
  if(all(fixed_effect != "")){
    line_noNA_Cofactor<-sapply(fixed_effect_split,function(x){which(!is.na(data[,x]))},simplify = F) %>% Reduce(intersect, .)
    NA_Cofactor<-nrow(data)-  length(line_noNA_Cofactor)
  }else{
    NA_Cofactor = 0
    line_noNA_Cofactor = 1:nrow(data)
  }
  
  if(NA_Cofactor != 0){
    warning(NA_Cofactor," observations containing missing values in the cofactors are elimnated.")
  }
  
  
  if(length(line_noNA_Y)!=nrow(data) | length(line_noNA_Cofactor)!=nrow(data) ){
    line_noNA<-intersect(line_noNA_Y,line_noNA_Cofactor)
    data<-data[line_noNA,]
    X<-X[line_noNA,]
  }
  
  #####################################################
  ##############        Response Y       ##############
  #####################################################
  
  #Extract Y 
  Y<-str_split(fixed,"[[~]]")[[1]][1] %>%
    str_replace_all(., " ", "")
  
  #check on Y
  if(Y==""){
    stop("No phenotype")
  }
  CheckColnames(data,Y)
  if(!class(data[,Y]) %in% c("numeric","interger","double")){
    stop("Y should be a numeric/interger/double")
  }
  
  Y<-data[,Y]
  
  #####################################################
  ##############      Fixed effects      ##############
  #####################################################
  
  #Saving Cofactors
  if(all(fixed_effect != "")){
    save_Cofactor<-data[,c(fixed_effect_split),drop=F]
  }
  
  #### classical effect without interaction
  
  #Which fixed effects are without interaction ?
  no_interaction<-fixed_effect[which(!fixed_effect %>%
                                       str_detect(., pattern = "\\:"))]
  
  #If mu is included in the modele, we add a new column mu full of 1 to data.
  if(any(no_interaction %in% "mu")){
    data$mu<-1
  }
  
  if(any(no_interaction %in% "" | length(no_interaction)==0)){ # if there are no fixed effects
    Liste_fixed<-NULL 
    
  }else{ #if there are fixed effects
    
    #Construction of indicator for each effect
    Liste_fixed<-lapply(no_interaction,function(x){CalcSimpleMatrixFixed(data,x)}) %>%
      map2(.x=.,.y=no_interaction, ~ {colnames(.x)<-paste0(.y,"_x_",colnames(.x));return(.x)})
    names(Liste_fixed)<-no_interaction
  }
  
  #### effect with interaction
  
  #Which fixed effects are with interactions ?
  interaction<-fixed_effect[which(fixed_effect %>%
                                    str_detect(., pattern = "\\:"))]
  
  
  if(length(interaction)==0){ # if there are no fixed effects with interaction
    Liste_MatrixInteraction<-NULL
    
  }else{ # if there are fixed effects with interaction
    
    #modifiction of interaction to access easily at each effect
    interaction<-lapply(interaction,function(x){FormatingInteraction(x)})
    
    #Construction of indicator for each effect
    Liste_MatrixInteraction<-lapply(interaction, function(x){CalcInteractionMatrixFixed(x,Liste_fixed)})
  }
  
  #### Compilation of fixed effect
  
  temp_names<-names(Liste_fixed)
  names(Liste_fixed)<-NULL
  
  if(is.null(Liste_fixed)){
    if(is.null(Liste_MatrixInteraction)){
      Cofactor<-NULL
    }else{
      Cofactor<-do.call("cbind",Liste_MatrixInteraction)
    }
  }else{
    if(is.null(Liste_MatrixInteraction)){
      Cofactor<-do.call("cbind",Liste_fixed)
    }else{
      Cofactor<-cbind(do.call("cbind",Liste_fixed),do.call("cbind",Liste_MatrixInteraction))
    }
  }
  names(Liste_fixed)<-temp_names
  
  #####################################################
  ##############      Error      ######################
  #####################################################
  
  #if thre is no error specify, we add a global error.
  if(!exists("error",.GlobalEnv)| str_remove_all(error,"~| ")==""){
    error<-"~ error"
  }
  
  #we add a column to data for error. Each observation will constitute a different factor. This allows us to consider error like any other random effect
  data$error<-factor(paste0("error",1:nrow(data)),levels=paste0("error",1:nrow(data)))
  
  #for the moment, K is defined as a diagonal matrix. -> modification to allow other matrix ?
  
  liste_error<-str_remove_all(error,"~| ") %>%
    str_split(.,"[[+]]") %>%
    .[[1]]
  
  Kerror<-diag(nrow(data))
  dimnames(Kerror)<-list(paste0("error",1:nrow(data)),paste0("error",1:nrow(data)))
  #because of the communication between function, we need to transfert Kerror in the global environnement.
  assign("Kerror", Kerror, envir = .GlobalEnv)
  
  error<-lapply(liste_error,function(x){
    if(!str_detect(x,",")){
      x<-str_replace_all(x,"error","error,Kerror")
    }
    if(!str_detect(x,"vm\\(")){
      x<-paste0("vm(",x,")")
    }
    return(x)
  }) %>%
    paste0(.,collapse=" + ")
  
  #####################################################
  ##############      Random effects     ##############
  #####################################################
  
  #special case : there are no random effect. But with still have an error that we will considered as a random effect.
  if(!exists("random")|str_remove_all(random," ")=="~"){
    random<-paste0("~",error)
  }else{
    random<-paste0(random," + ",error)
  }
  
  #### effect with no interaction
  
  #extract random effect names
  random_effect<-str_split(random,"[[~]]")[[1]][2] %>%
    str_replace_all(., " ", "") %>%
    str_split(.,"[[+]]")%>%
    .[[1]] %>%
    str_remove_all(.,"vm\\(|\\)")
  
  
  
  #Which random effects are without interactions ?
  no_interaction<-random_effect[which(!random_effect %>%
                                        str_detect(., pattern = "\\:"))] %>%
    str_split(.,",")
  
  
  if(length(no_interaction)==0){# if there are no random effects without interaction
    Liste_random<-NULL
  }else{# if there are random effects without interaction
    
    #Construction of design matrix and variance covariance matrix for each effect
    Liste_random<-lapply(no_interaction,function(x){CalcSimpleMatrixRandom(data,x)})
    names(Liste_random)<-lapply(no_interaction,function(x){paste0(x,collapse=",")}) %>%
      paste0("(",.,")")
  }
  
  
  #### effect with interaction
  
  #Which random effects are with interactions ?
  interaction<-random_effect[which(random_effect %>%
                                     str_detect(., pattern = "\\:"))] %>%
    str_split(.,":")
  
  if(length(interaction)==0){# if there are no random effects with interaction
    Liste_random_interaction<-NULL
  }else{# if there are random effects with interaction
    
    #Construction of design matrix and variance covariance matrix for each effect
    Liste_random_interaction<-lapply(interaction,function(x){CalcInteractionMatrixRandom(data,x,Liste_random,Liste_fixed)})
    
  }
  
  #####################################################
  ##############   Compilation   ######################
  #####################################################
  
  #compilation of all design matrix (ZList) and variance-covariance matrix (VList)
  VList<-compilation_random(Liste_random,Liste_random_interaction,interaction,Type="V")
  ZList<-compilation_random(Liste_random,Liste_random_interaction,interaction,Type="Z")
  
  detect_error<-str_detect(names(VList),"error")
  VList<-VList[c(which(!detect_error),which(detect_error))]
  ZList<-ZList[c(which(!detect_error),which(detect_error))]
  
  
  #Cleaning of the global environnement
  
  if(exists("Kerror")){
    rm(Kerror,envir=.GlobalEnv)
  }
  
  #rm(data, envir = .GlobalEnv)
  
  #Formating of Formula
  Formula<-str_split(fixed," |~|[[+]]") %>% unlist() %>% .[-which(. == "")] %>% .[-1]
  
  if(length(Formula)==0){
    Formula<-NULL
    save_Cofactor<-NULL
  }else{
    Formula <- Formula %>%.[which(.!="mu")] %>% paste0(.,collapse="+") %>% paste0("~",.) %>% as.formula()
  }
  
  
  
  return(list(Y=Y,Cofactor=save_Cofactor,X=X,Formula=Formula,ZList=ZList,VList=VList))
  
}








###########################################################################################################
#####################               MultipleEasyUseMM4LMM             #####################################
###########################################################################################################

FormuleEnsemble<-function(x,trait,element){
  
  if(element == "fixed"){
    sapply(x,function(y){str_remove(List_modele[[y]][element],paste0(trait," ~ "))}) %>%
      paste(.,collapse = " + ") %>%
      str_split(.," [[+]] ") %>% 
      .[[1]] %>%
      str_remove_all(.," ") %>%
      unique(.) %>%
      paste(.,collapse = " + ") %>%
      paste0(trait," ~ ",.)
  }else{
    
    sapply(x,function(y){if(element %in% names(List_modele[[y]])){str_remove(List_modele[[y]][element],"~ ") %>% str_remove(.,"~")}else{""}}) %>%
      .[which(. != "")] %>%
      paste(.,collapse = " + ") %>%
      str_split(.," [[+]] ") %>% 
      .[[1]] %>%
      str_remove_all(.," ") %>%
      unique(.) %>%
      paste(.,collapse = " + ") %>%
      paste0("~ ",.)
    
    
  }
}

SelectionMatrixByModelByELement<-function(x,Matrice_ensemble,trait,element){
  if(element == "fixed"){
    Matrix_type = "Cofactor"
    Name_matrix_type<-names(Matrice_ensemble[[trait]][[Matrix_type]])
    
    Selected_matrix<-x[element] %>% str_split(.," ~ ") %>% .[[1]] %>% .[2] %>% str_split(., " [[+]] ") %>% .[[1]] %>% 
      str_remove_all(.," |at\\(|vm\\(|\\)") %>%
      sapply(.,function(y){str_split(y,":")}) %>%
      lapply(.,function(y){y<-str_replace(y,",","_x_");
      if(length(y) == 1){
        which(str_detect( Name_matrix_type,y) 
              & !str_detect(Name_matrix_type,":")) %>%
          Name_matrix_type[.] %>% unique(.)
      }else{
        which(str_detect( Name_matrix_type,y[1]) & str_detect( Name_matrix_type,y[2])
              & str_detect(Name_matrix_type,":")) %>%
          Name_matrix_type[.] %>% unique(.)
      }}) %>% unlist(.) %>% unique(.)
    return(Matrice_ensemble[[trait]][["Cofactor"]][,Selected_matrix])
    
  }else{
    Matrix_type = "ZList"
    Name_matrix_type<-names(Matrice_ensemble[[trait]][[Matrix_type]])
    
    Selected_matrix<-x[element] %>% str_split(., " [[+]] ") %>% .[[1]] %>% 
      str_remove_all(.,"~| |vm\\(|\\)") %>%
      sapply(.,function(y){str_split(y,":")}) %>%
      lapply(.,function(y){
        
        if(length(y) == 1){#no interaction
          if(!str_detect(y[1],"error,")){
            y[1]<-str_replace(y[1],"error",fixed("error,Kerror"))
          }
          y[1]<-paste0("\\(",y[1],"\\)")
          which(str_detect( Name_matrix_type,y) & !str_detect(Name_matrix_type,":")) %>% Name_matrix_type[.]
          
        }else{#interaction
          
          if(!str_detect(y[2],"error,")){#error without kinship
            y[2]<-str_replace(y[2],"error",fixed("error,Kerror"))
          }
          y[1]<-str_replace(y[1],",","_x_")
          y[2]<-paste0("(",y[2],")")
          if(str_detect(y[1],"_x_")){#Specific levels
            y[1]<-str_remove(y[1],"at\\(")
            which(str_detect( Name_matrix_type,y[1]) & str_detect( Name_matrix_type,y[2])
                  & str_detect(Name_matrix_type,":")) %>% 
              Name_matrix_type[.] %>% unique(.)
            
          }else{#Global
            if(str_detect(y[1],"at\\(")){#one matrix by level
              y[1]<-str_remove(y[1],"at\\(")
              which(str_detect( Name_matrix_type,y[1]) & str_detect( Name_matrix_type,fixed(y[2])) 
                    & str_detect(Name_matrix_type,":") & str_detect(Name_matrix_type,"_x_")) %>% 
                Name_matrix_type[.]  %>% unique(.)
              
            }else{#one matrix
              which(str_detect( Name_matrix_type,y[1]) & str_detect( Name_matrix_type,fixed(y[2])) 
                    & str_detect(Name_matrix_type,":") & !str_detect(Name_matrix_type,"_x_")) %>% 
                Name_matrix_type[.]  %>% unique(.)
            }
          }
        }}) %>% unlist(.) %>% unique(.)
    
    return(list(ZList=Matrice_ensemble[[trait]][["ZList"]][Selected_matrix],
                VList=Matrice_ensemble[[trait]][["VList"]][Selected_matrix]))
  }
  
}

SelectionMatrixByModel<-function(x,Matrice_ensemble){
  Matrix_selected<-list()
  trait<-str_split(x$fixed, "[[~]]") %>% .[[1]] %>% .[1] %>% str_remove_all(.," ")
  
  #Y response
  Matrix_selected[["Y"]]<-Matrice_ensemble[[trait]][["Y"]]
  
  #fixed part
  Matrix_selected[["Cofactor"]]<-SelectionMatrixByModelByELement(x,Matrice_ensemble,trait,element="fixed")
  
  #ZList, VList
  Random_part<-SelectionMatrixByModelByELement(x,Matrice_ensemble,trait,element="random")
  if("error" %in% names(x)){
    Error_part<-SelectionMatrixByModelByELement(x,Matrice_ensemble,trait,element="error")
  }else{
    Error_part<-list()
  }
  
  Matrix_selected[["ZList"]]<-c(Random_part$ZList,Error_part$ZList)
  Matrix_selected[["VList"]]<-c(Random_part$ZList,Error_part$VList)
  
  return(Matrix_selected)
}


MultipleEasyUseMM4LMM<-function(data=data,List_modele=List_modele){
  
  #Check element names
  if(!all(sapply(names(List_modele),function(x){all(names(List_modele[[x]]) %in% c("fixed","random","error"))}))){
    stop("For each model, element names should be 'fixed','random' or 'error'")
  }
  
  if(!all(sapply(names(List_modele),function(x){'fixed' %in% names(List_modele[[x]])}))){
    stop("For each model, you should indicate at least 'fixed' if there is any random effects use '~ '")
  }
  
  List_modele<-lapply(List_modele,function(y){if(!"error" %in% names(y)){y[["error"]]<-"~ vm(error)"};y})
  List_modele<-lapply(List_modele,function(y){if(y[["error"]] == "~ "){y[["error"]]<-"~ vm(error)"};y})
  List_modele<-lapply(List_modele,function(y){if(!"random" %in% names(y)){y[["random"]]<-"~ "};y})
  
  lapply(List_modele, function(x){CheckNamesModele(data,fixed=x[["fixed"]],random=x[["random"]],error=x[["error"]])})
  
  #Check response
  if(!all(unlist(lapply(List_modele,function(x){str_split(x[["fixed"]],"~")[[1]][1] %>% str_split(.," ") %>% .[[1]] %>% length(.) == 2})))){
    stop("Please, Check response part in fixed part")
  }
  
  #Modele by trait
  traits<-lapply(List_modele,function(x){str_split(x[["fixed"]],"~")[[1]][1] %>% str_split(.," ") %>% .[[1]] %>% .[1]})
  Modele_by_trait<-sapply(unique(unlist(traits)),function(x){which(traits == x) %>% traits[.] %>% names(.)},simplify=F)
  
  
  #Global Model
  fixed_ensemble<-imap(Modele_by_trait,function(.x,.y){FormuleEnsemble(x=.x,trait=.y,element='fixed')})
  random_ensemble<-imap(Modele_by_trait,function(.x,.y){FormuleEnsemble(x=.x,trait=.y,element='random')})
  error_ensemble<-imap(Modele_by_trait,function(.x,.y){FormuleEnsemble(x=.x,trait=.y,element='error')})
  
  Matrice_ensemble<-sapply(names(Modele_by_trait),function(x){EasyUseMM4LMM(data=pheno,fixed= fixed_ensemble[[x]],random=random_ensemble[[x]],error=error_ensemble[[x]])},simplify=F)
  
  #Matrices for each modele
  Matrix_by_model<-lapply(List_modele, function(x){SelectionMatrixByModel(x,Matrice_ensemble)})
  
  return(Matrix_by_model)
}



