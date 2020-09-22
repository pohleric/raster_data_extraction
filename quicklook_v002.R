
#---------------------------------------------- FUNCTIONS ----------------------------------------------
prep <- function(fname,adj_fac=1, offset = 0){
  # fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_GPM-3IMERG_mon/precipitation_1.txt'
  lmfun <- function(t,x){
    LM <- lm(x~t)
    anoms <- LM$residuals
    return(anoms)
  }
  print(paste0('Reading ',fname))
  a <- read.table(fname, sep=',', header = T)
  colnames(a)[1] <- 'time'
  a$time <- as.POSIXct(a$time, tz="UTC")+6*60*60 # Kyrgyz time = UTC+6
  # a$time <- as.POSIXct(a$time)
  # n_col = NCOL(a) - 1
  names_c <- colnames(a)
  names_c <- names_c[names_c!='time']
  a[names_c] <- (a[names_c] * adj_fac) - offset
  
  a_cop = a
  ind_0 = is.na(a_cop[names_c])
  ind_val = !is.na(a_cop[names_c])
  
  a_cop[is.na(a_cop)] <- 0
  
  # go through each column (except time) and calculate the cumulative sum and the detrended cumulative sum
  print(paste0('Calculating  cumulative sum and anomalies (with respect to the cumulative sums)'))
  tmp_cs <- apply(a_cop[names_c], 2, cumsum)
  tmp_cs[ind_0] <- NA
  # head(tmp_cs)
  # colnames(tmp_cs)
  tmp_cs_dt <- apply(tmp_cs, 2, function(x)lmfun(a$time,x))
  
  # write the cs and cs_dt back into dataframe
  names_c_cs <- paste0(names_c,'_cs')
  names_c_cs_dt <- paste0(names_c_cs,'_dt')
  a[names_c_cs] <- tmp_cs
  a[names_c_cs_dt] <- NA
  ind_val_multi = (rowMeans(ind_val)==1)
  a[ind_val_multi,names_c_cs_dt] <- tmp_cs_dt
  
  # a[as.logical(ind_0),names_c_cs] <- NA
  # a[as.logical(ind_0),names_c_cs_dt] <- NA
  
  
  print('Done.')
  return(a)
}

# time series format?
zoo_plot <- function(x,what='cs_dt', average =FALSE){
  par(mfrow=c(1,1))
  # x <- d
  require('zoo')
  # x <- d
  data_zoo <- zoo(x[,2:NCOL(x)],order.by = x$time)
  # cs_dt plot
  if(what=='cs_dt'){
    ind_cs_dt <- grep('cs_dt$', names(data_zoo))
  }
  if(what=='cs'){
    ind_cs_dt <- grep('cs$', names(data_zoo))
  }
  
  if(what=='ts'){
    # ind_cs_dt <- substr(names(data_zoo))
    nch <- nchar(names(data_zoo))
    lastc <- substr(names(data_zoo),start = nch, stop = nch)
    lastc_num <- as.numeric(lastc)
    ind_cs_dt <- which(!is.na(lastc_num))
  }
  # if too many points - calc mean and quantile bounds
  if(average == TRUE){
    data__mean <- apply(data_zoo[,ind_cs_dt],1,mean)
    data__qnt25 <- apply(data_zoo[,ind_cs_dt],1,quantile,probs=c(0.25))
    data__qnt75 <- apply(data_zoo[,ind_cs_dt],1,quantile,probs=c(0.75))
    
    data_zoo_mean <- zoo(data__mean, order.by = as.POSIXct(names(data__mean), tz='UTC'))
    data_zoo_q25 <- zoo(data__qnt25, order.by = as.POSIXct(names(data__qnt25), tz='UTC'))
    data_zoo_q75 <- zoo(data__qnt75, order.by = as.POSIXct(names(data__qnt75), tz='UTC'))
    # t(data__qnt)
    
    plot(data_zoo_mean, xlab = names(x)[1], ylab = paste0('average_',what))
    lines(data_zoo_q25,lty=3)
    lines(data_zoo_q75,lty=3)
    abline(h=0, col='red', lty=2)
  }else{
    n <- length(ind_cs_dt)
    rows <- ceiling(sqrt(n))
    cols <- floor(sqrt(n))
    par(mfrow=c(rows, cols))
    for(ind_i in ind_cs_dt){
      print( names(x)[ind_i])
      plot(data_zoo[,ind_i], xlab = names(x)[1], ylab = names(data_zoo)[ind_i])
      abline(h=0, col='red', lty=2)
    }
  }
}

#---------------------------------------------- DATA ----------------------------------------------


fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_GPM-3IMERG_mon/precipitation_1.txt'
d <- prep(fname)
zoo_plot(d, what='ts')
zoo_plot(d, what = 'cs')
zoo_plot(d, what = 'cs_dt', average=T)


# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_GPCC/precip_1.txt'
# d <- prep(fname)
# zoo_plot(d, what='ts')
# zoo_plot(d, what = 'cs')
# zoo_plot(d, what = 'cs_dt',average = T)


# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_APHR_V1101/precip_1.txt'
# d <- prep(fname)
# # zoo_plot(d, what='ts')
# # zoo_plot(d, what = 'cs')
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_APHR_V1101_EXR1/precip_1.txt'
# d <- prep(fname)
# zoo_plot(d,what = 'cs_dt')
# 
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/tp_1.txt'
# d <- prep(fname, adj_fac=1000)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_APHR_V1808/tave_1.txt'
# d <- prep(fname, offset = -273.15)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_APHR_V1808R1/precip_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'ts')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_chelsacruts_mon/prec_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt', average=TRUE)
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_chelsats_mon/prec_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt', average=TRUE)
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_chelsats_mon/tmean_1.txt'
# d <- prep(fname, adj_fac=1/10)
# zoo_plot(d, what = 'ts', average=TRUE)
# # temperature is too high by a factor 10
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/d2m_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'ts')
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/sp_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'ts')
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/ssr_1.txt'
# d <- prep(fname)
# # zoo_plot(d, what = 'ts')
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/t2m_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/tcc_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/tp_1.txt'
# d <- prep(fname, adj_fac = 1000)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/u10_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/v10_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR10_V2_day/lwdown_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR10_V2_day/lwup_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR10_V2_day/prcp_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR10_V2_day/q2_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR10_V2_day/swdown_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR10_V2_day/swup_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR10_V2_day/t2_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR10_V2_day/ws10_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_day/lwdown_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_day/lwup_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_day/prcp_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_day/psfc_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_day/q2_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_day/swdown_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_day/swup_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')
# 
# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_day/t2_1.txt'
# d <- prep(fname)
# zoo_plot(d, what = 'cs_dt')

# fname <- '/Users/pohle/PycharmProjects/data-download/marlene/output_CRU-TS/pre_1.txt'
# d <- prep(fname)
# # zoo_plot(d, what = 'ts')
# zoo_plot(d, what = 'cs_dt')
