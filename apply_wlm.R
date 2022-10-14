#!/usr/binenv Rscript
# Scipt to load R functions to apply WLM in pairs or trios with arbitraty sample overlap
# For errors or bugs contact Robin N Beaumont r.beaumont@exeter.ac.uk
print("")

apply_wlm_trio<-function(data,cm,cp,mp){
  # trio
  data$c_wlm_beta<-2*data$c_beta-data$m_beta-data$p_beta
  data$c_wlm_se<-sqrt(4*data$c_se**2+data$m_se**2+data$p_se**2-4*cm*sqrt((data$c_se**2)*(data$m_se**2))-4*cp*sqrt((data$c_se**2)*(data$p_se**2))+2*mp*sqrt((data$m_se**2)*(data$p_se**2)))
  data$c_wlm_p<-2*pnorm(-abs(data$c_wlm_trio_beta/data$c_wlm_trio_se))
  data$m_wlm_beta<-(3*data$m_beta-2*data$c_beta+data$p_beta)/2
  data$m_wlm_se<-sqrt((9/4)*data$m_se**2+data$c_se**2+(1/4)*data$p_se**2-3*cm*sqrt((data$c_se**2)*(data$m_se**2))-cp*sqrt((data$c_se**2)*(data$p_se**2))+1.5*mp*sqrt((data$m_se**2)*(data$p_se**2)))
  data$m_wlm_p<-2*pnorm(-abs(data$m_wlm_trio_beta/data$m_wlm_trio_se))
  data$p_wlm_beta<-(3*data$p_beta-2*data$c_beta+data$m_beta)/2
  data$p_wlm_se<-sqrt((9/4)*data$p_se**2+data$c_se**2+(1/4)*data$m_se**2-3*cp*sqrt((data$c_se**2)*(data$p_se**2))-cm*sqrt((data$c_se**2)*(data$m_se**2))+1.5*mp*sqrt((data$m_se**2)*(data$p_se**2)))
  data$p_wlm_p<-2*pnorm(-abs(data$p_wlm_trio_beta/data$p_wlm_trio_se))
  # Return data frame
  data
}

apply_wlm_pair<-function(data,cm){
  # pair
  data$c_wlm_beta<-(4*data$c_beta-2*data$m_beta)/3
  data$c_wlm_se<-sqrt((16*data$c_se**2+4*data$m_se**2-16*cm*sqrt((data$c_se**2)*(data$m_se**2)))/9)
  data$c_wlm_p<-2*pnorm(-abs(data$c_wlm_beta/data$c_wlm_se))
  data$m_wlm_beta<-(4*data$m_beta-2*data$c_beta)/3
  data$m_wlm_se<-sqrt((16*data$m_se**2+4*data$c_se**2-16*cm*sqrt((data$c_se**2)*(data$m_se**2)))/9)
  data$m_wlm_p<-2*pnorm(-abs(data$m_wlm_beta/data$m_wlm_se))
  # Return data frame
  data
}

# data<-apply_wlm_trio(data,cm,cp,mp)
