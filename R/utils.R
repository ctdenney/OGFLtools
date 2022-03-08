.onAttach <- function(libname, pkgname) {
  packageStartupMessage("><(((O> ><(((O> ><(((O> Welcome to the OGFL utility package <O)))>< <O)))>< <O)))>< ")

}

sn = function(x) {
  length(unique(x))
}

ec_to_sal = function(temp, cond){
  if (any(temp > 35, na.rm = T)) {
    warning('Temperature is high, ensure that units are in degrees C', call. = F, immediate. = T)
  }
  if(length(temp == 1)) {
    temp = rep(temp, length.out = length(cond))
  } else if ( length(temp == length(cond))) {
    temp = temp
  } else {
    warning('Temperature must be of length 1 or match number of conductivity readings')
  }
  ref_cond = 42914
  cond_rat = cond/ref_cond
  rt = 0.6766097 + (0.0200564*temp) + (0.0001104259*(temp^2)) + ((-6.9698*10^-7)*(temp^3)) + ((1.0031*10^-9)*temp^4)
  Rt = cond_rat/rt
  dS = ((temp-15)/(1+0.0162*(temp-15)))*(0.0005+(-0.0056)*(Rt^0.5)+(-0.0066)*Rt+(-0.0375)*(Rt^1.5)+(0.0636)*(Rt^2)+(-0.0144)*(Rt^2.5))
  sal = c()
  sal[is.na(cond)] = NA
  sal[cond >3000 & !is.na(cond)] = 0.008 + (-0.1692)*(Rt[cond>3000 & !is.na(cond)]^0.5) + 25.3851*Rt[cond > 3000 & !is.na(cond)]+14.0941*(Rt[cond > 3000 & !is.na(cond)]^1.5)+(-7.0261)*(Rt[cond > 3000 & !is.na(cond)]^2)+2.7081*(Rt[cond > 3000 & !is.na(cond)]^2.5)+dS[cond > 3000 & !is.na(cond)]
  sal[cond <=3000 & !is.na(cond)] = (0.008+ (-0.1692)*(Rt[cond <= 3000 & !is.na(cond)]^0.5)+25.3851*Rt[cond <= 3000 & !is.na(cond)]+14.0941*(Rt[cond <= 3000 & !is.na(cond)]^1.5)+(-7.0261)*(Rt[cond <= 3000 & !is.na(cond)]^2)+2.7081*(Rt[cond <= 3000 & !is.na(cond)]^2.5)+dS[cond <= 3000 & !is.na(cond)]) -
    (0.008/(1+(1.5*(400*Rt[cond <= 3000 & !is.na(cond)]))+((400*Rt[cond <= 3000 & !is.na(cond)])^2))-(0.0005*(temp[cond <= 3000 & !is.na(cond)]-15)/(1+0.0162*(temp[cond <= 3000 & !is.na(cond)]-15)))/(1+((100*Rt[cond <= 3000 & !is.na(cond)])^0.5)+((100*Rt[cond <= 3000 & !is.na(cond)])^1.5)))
  return(sal)
}

sal_to_spc = function(sal){
  if (any(sal > 40, na.rm = T)){
    warning("Salinity is high, ensure that data is correct")
  }
  J1 = -16.072
  J2 = 4.1495
  J3 = -0.5345
  J4 = 0.0261
  spc = (sal/35)*(53087) + sal*(sal-35) *
    (J1 + (J2 * (sal^(0.5))) + (J3*sal) + (J4*(sal^(1.5))))
  return(spc)
}


# sr_to_sal = function(sr, srfw = 0.705264, srmar = 0.70918,confw = 74.6, conmar = 6819,salfw = 0.1,salmar = 31.8, sallim_high, fill_high = "NA", sallim_low, fill_low = "NA",suppress.warnings = F){
#   if (any(sr < min(srfw, srmar)| sr > max(srfw, srmar), na.rm = T) & suppress.warnings == F) {
#     warning('Some of your measured strontium ratio values are outside the bounds of your two endmembers, make sure that srfw and srmar are set correctly',
#             call. = F, immediate. = T)
#   }
#   if (is.na(sallim_high) | is.na(sallim_low)) stop('You have not set one or both of the salinity limit values(sallim_high and sallim_low arguments)')
#     sal = (((salfw*srmar*conmar) - (salfw*sr*conmar) - (salmar*srmar*conmar) + (salmar*sr*conmar))/
#              ((sr*confw) - (sr*conmar) - (srfw*confw) + (srmar*conmar))) + salmar
#     if((fill_high != 'NA' & is.numeric(fill_high) == F) | (fill_low != 'NA' & is.numeric(fill_low) == F)) stop('you have not set an appropriate fill value. Fill must be either \'NA\' or a numeric value, for both fill_low and fill_high.')
#     else if(srmar > srfw){
#       if (fill_high == 'NA') {
#         sal[sal > sallim_high | sr > srmar] = NA
#       } else if (is.numeric(fill_high)) {
#         sal [sal > sallim_high |sr > srmar] = fill_high
#       }
#
#       if(fill_low == 'NA') {
#         sal[sal < sallim_low | sr < srfw] = NA
#       } else if (is.numeric(fill_low)) {
#         sal [sal < sallim_low | sr < srfw] = fill_low
#       }
#     } else if (srmar < srfw) {
#       if (fill_high == 'NA') {
#         sal[sal > sallim_high | sr < srmar] = NA
#       } else if (is.numeric(fill_high)) {
#         sal [sal > sallim_high |sr < srmar] = fill_high
#       }
#
#       if(fill_low == 'NA') {
#         sal[sal < sallim_low | sr > srfw] = NA
#       } else if (is.numeric(fill_low)) {
#         sal [sal < sallim_low | sr > srfw] = fill_low
#       }
#     }
#     return(sal)
# }

sr_to_sal = function(sr, sr1 = 0.70576, sr2 = 0.70918,con1 = 95, con2 = 6819,sal1 = 0.11,sal2 = 31.8, sallim_high = max(sal1, sal2), fill_high = "NA", sallim_low = min(sal1,sal2), fill_low = "NA",suppress.warnings = F){

  if (sallim_low < min(sal1, sal2) | sallim_high > max(sal1,sal2)) {
    stop('One or both of your sallims are outside the bounds of your endmembers')
  }
  if (is.na(sallim_high) | is.na(sallim_low)) stop('You have not set one or both of the salinity limit values(sallim_high and sallim_low arguments)')
  if (any(sr < min(sr1, sr2)| sr > max(sr1, sr2), na.rm = T) & suppress.warnings == F) {
    warning('Some of your measured strontium isotope ratio values are outside the bounds of your two endmembers, make sure that sr1 and sr2 are set correctly',
            call. = F, immediate. = T)
  }
  sal = (((sal1*sr2*con2) - (sal1*sr*con2) - (sal2*sr2*con2) + (sal2*sr*con2))/
           ((sr*con1) - (sr*con2) - (sr1*con1) + (sr2*con2))) + sal2
  if((fill_high != 'NA' & is.numeric(fill_high) == F) | (fill_low != 'NA' & is.numeric(fill_low) == F)) stop('you have not set an appropriate fill value. Fill must be either \'NA\' or a numeric value, for both fill_low and fill_high.')
  else if(sr2 > sr1){
    if (fill_high == 'NA') {
      sal[sal > sallim_high | sr > sr2] = NA
    } else if (is.numeric(fill_high)) {
      sal [sal > sallim_high |sr > sr2] = fill_high
    }

    if(fill_low == 'NA') {
      sal[sal < sallim_low | sr < sr1] = NA
    } else if (is.numeric(fill_low)) {
      sal [sal < sallim_low | sr < sr1] = fill_low
    }
  } else if (sr2 < sr1) {
    if (fill_high == 'NA') {
      sal[sal > sallim_high | sr < sr2] = NA
    } else if (is.numeric(fill_high)) {
      sal [sal > sallim_high | sr < sr2] = fill_high
    }

    if(fill_low == 'NA') {
      sal[sal < sallim_low | sr > sr1] = NA
    } else if (is.numeric(fill_low)) {
      sal [sal < sallim_low | sr > sr1] = fill_low
    }
  }
  return(sal)
}

sal_to_sr = function(sal, sr1 = 0.70576, sr2 = 0.70918, con1 = 95, con2 = 6819, sal1 = 0.11, sal2 = 31.8) {
  if(any(sal < min(sal1, sal2) | sal > max(sal1, sal2), na.rm = T)) {
    warning('Some of your measured salinity values are outside the bounds of your two endmembers, make sure that sal1 and sal2 are set correctly',
            call. = F, immediate. = T)
  }
    sr = ((((sr1*con1*sal)-(sr1*con1*sal2))/(sal1 - sal2))+(sr2*con2)-(((sr2*con2*sal)-(sr2*con2*sal2))/(sal1 - sal2)))/
      ((((con1*sal)-(con1*sal2))/(sal1 - sal2))+(con2)-(((con2*sal)-(con2*sal2))/(sal1-sal2)))
    return(sr)

}

o_to_sal = function(d18O, source = 'ingram') {
  if(source == 'ingram') {
    sal = (d18O + 10.9)/0.32
  } else if(source == 'mclg') {
    sal = (d18O + 10.17)/0.29
  } else{
    warning('Non-supported source')
  }
  return(sal)
}

vsmow_to_vpdb = function(x) {
  return((0.97001*x)-29.99)
}

vpdb_to_vsmow = function(x) {
  return((1.03091*x)+30.91)
}

bim = function(fl, hl, gt) {
  fl + ((gt - max(gt))*(fl - hl))/(max(gt) - min(gt))
}

membermix = function(sr, conc, sal, mix) {
  if(sum(mix) != 1) {
    warning('Your mixture does not sum to 100%',
            call. = F, immediate. = T)
  }
  srmix = sum(sr*conc*mix)/sum(conc*mix)
  srconc = sum(conc*mix)
  salmix = sum(sal*mix)
  return(list(sr = srmix, conc = srconc, sal = salmix))
}

##all these values come from models based on the lfs_length_data.rda file included with this package
##the models are calculated in the lengthconversion script
length_con = function(from, to = 'flf',length) {
  if(from == 'slf' & to == 'flf') {
    calclength =  1.046694 *(length)+ 2.741875
  } else if (from == 'slf' & to == 'tlf') {
    calclength =  1.164434 *length+ 1.095035
  } else if (from == 'tlf' & to == 'flf') {
    calclength =  0.8959148 *length+ 2.030626
  } else if (from == 'tlf' & to == 'slf') {
    calclength =  0.8530776 *length+ -0.4158636
  } else if (from == 'flf' & to == 'slf') {
    calclength =  0.9495577 *length+ -2.127712
  } else if (from == 'flf' & to == 'tlf') {
    calclength =  1.109416 *length+ -1.696259
  } else if (from == 'sle' & to == 'flf') {
    ## used sle-fle-flf instead of sle-slf-flf because the sle-slf model is considerably worse than the sle-fle model
    intercalc =  1.092 *length+ 0.8059086 #sle to fle
    calclength =  1.017514 *intercalc+ 1.492639 #fle to flf
  } else if (from == 'sle' & to == 'tlf') {
    intercalc =  1.023139 *length+ 3.570214
    calclength =  1.164434 *intercalc+ 1.095035
  } else if (from == 'sle' & to == 'slf') {
    calclength =  1.023139 *length+ 3.570214
  } else if (from == 'fle' & to == 'flf') {
    calclength =  1.017514 *length+ 1.492639
  } else if (from == 'fle' & to == 'slf') {
    intercalc =  1.017514 *length+ 1.492639 # fle to flf
    calclength =  0.9495577 *intercalc+ -2.127712 #flf to slf
  } else if (from == 'fle' & to == 'tlf') {
    intercalc =  1.017514 *length+ 1.492639 #fle to flf
    calclength =  1.109416 *intercalc+ -1.696259 #flf to tlf
  } else if (from == 'flfrz' & to == 'flf') {
    calclength =  0.9983984 *length+ 2.084482
  } else if (from == 'slfrz' & to == 'flf') {
    calclength =  1.052409 *length+ 4.533137
  } else if (from == 'tlfrz' & to == 'flf') {
    calclength =  0.9055277 *length+ 2.737323
  }else {
    stop('Unsupported conversion', call. = F)
  }
  if(any(length < 20 | length > 120, na.rm = T)) {
    warning("Some lengths outside of established conversion bounds", call. = F)
  }
  return(round(calclength,1))
}


colmatch <- function(x,y, match = F, join = F){
  if (match == T & join == F) {
    if (all(colnames(x) %in% colnames(y)) & all(colnames(y) %in% colnames(x))) {
      return('All column names match between both objects')
    } else if (any(colnames(x) %in% colnames(y))) {
      cat("\n\n",'columns ')
      cat(colnames(x)[colnames(x) %in% colnames(y)], sep = " and ")
      cat( " are found in both objects\n")
    } else if(!any(colnames(x) %in% colnames(y))) {
      cat("\n\n",'No columns ')
      cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
      cat( " are matching between both objects\n")
    }
  } else if(match == F & join == F) {
    if (all(colnames(x) %in% colnames(y)) & all(colnames(y) %in% colnames(x))) {
      return('All column names match between both objects')
    } else if (!all(colnames(x) %in% colnames(y))
               & all(colnames(y) %in% colnames(x))) {
      if(length(colnames(x)[!colnames(x) %in% colnames(y)]) > 1) {
        cat("\n\n",'columns ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " are missing from second object\n")
      } else {
        cat("\n\n",'column ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " is missing from second object\n")
      }
    } else if (all(colnames(x) %in% colnames(y)) & !all(colnames(y) %in% colnames(x))){
      if(length(colnames(y)[!colnames(y) %in% colnames(x)]) > 1) {
        cat("\n\n",'columns ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat( " are missing from first object\n")
      } else {
        cat("\n\n",'column ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat( " is missing from first object\n")
      }
    } else if (!all(colnames(x) %in% colnames(y)) & !all(colnames(y) %in% colnames(x))){
      if(length(colnames(y)[!colnames(y) %in% colnames(x)]) > 1 & length(colnames(x)[!colnames(x) %in% colnames(y)]) > 1) {
        cat("\n\n",'columns ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat(' are missing from first object \n and')
        cat('\n','columns ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " are missing from second object\n")
      } else if(length(colnames(y)[!colnames(y) %in% colnames(x)]) > 1 & length(colnames(x)[!colnames(x) %in% colnames(y)]) == 1){
        cat("\n\n",'columns ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat(' are missing from first object \n and')
        cat('\n','column ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " is missing from second object\n")
      } else if(length(colnames(y)[!colnames(y) %in% colnames(x)]) == 1 & length(colnames(x)[!colnames(x) %in% colnames(y)]) > 1){
        cat("\n\n",'column ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat(' is missing from first object \n and')
        cat('\n','columns ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " are missing from second object\n")
      } else {
        cat("\n\n",'column ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat(' is missing from first object \n and')
        cat('\n','column ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " is missing from second object\n")
      }
    }
  } else if (match == T & join == T) {
    colnames(x[colnames(x) %in% colnames(y)])
    } else {
    cat('Imporper selection of match argument, select either "TRUE/T" or "FALSE/F. If join is set to true, match must also be T)')
  }

}

ec_to_spc = function(temp, ec) {
  spc = ec/(1 + (0.0191 * (temp - 25)))
  return(spc)
}

spc_to_ec = function(temp,spc) {
  ec = (1 + (0.0191 * (temp - 25)))*spc
  return(ec)
}
