library(chunked)
library(tabplot)
options(fftempdir = "~/rfiles/tmp")

indir <- "/media/disk1/Eucalyptus_temp/singleSp_hwe/hwe_merge"
in_file <- file.path(indir,"combined.hwe")
#in_file <- file.path(getwd(),"tmp2")
file_chunked<-
  read_chunkwise(in_file,format="table")
head(file_chunked,n=5L) # %>% View()
x=file_chunked %>% select(chromo) %>% head(n=5L)
table(x)


suffix=".hwe.Major"
species=file_chunked %>% head() %>% names() %>% grep(suffix,.,value=T) %>% gsub(suffix,"",.)


Ho_calc <- function(freq,Fis){
  2*freq*(1-freq)-2*freq*(1-freq)*Fis
}

speciesHo <- function(df,spChar){
  require("dplyr")
  require("lazyeval")
  freqName=paste(spChar,".hwe.Freq",sep="")
  FisName=paste(spChar,".hwe.F",sep="")
  HoName=paste(spChar,".hwe.Ho",sep="")
  mutate_call = lazyeval::interp(~ Ho_calc(a,b), a = as.name(freqName), b = as.name(FisName))
  df %>% mutate_(.dots = setNames(list(mutate_call), HoName))
}
#head(speciesHo(speciesHo,"caleyi_common"))

  maxHo=0.8

# file_chunked %>%
#   #  select(one_of(paste(species,suffix,sep=""))) %>%
#   apply{
#
#   }
# head()

myframe <- file_chunked %>%
  select(snpID, chromo, position,
         albens_common.hwe.Freq,albens_common.hwe.F,albens_common.hwe.p.value,
         sideroxylon_common.hwe.Freq,sideroxylon_common.hwe.F,sideroxylon_common.hwe.p.value) %>%
  # select(snpID, chromo, position, starts_with("albens"), starts_with("sideroxylonal")) %>%
  # head() %>%
        mutate(albens_common.hwe.Ho = Ho_calc(albens_common.hwe.Freq,albens_common.hwe.F),
         sideroxylon_common.hwe.Ho = Ho_calc(sideroxylon_common.hwe.Freq,sideroxylon_common.hwe.F),
         albens_common_para = albens_common.hwe.Ho > 0.7,
         sideroxylon_common_para = sideroxylon_common.hwe.Ho > 0.7,
         albens_common_hwe01 = albens_common.hwe.p.value < 0.01,
         sideroxylon_common_hwe01 = sideroxylon_common.hwe.p.value < 0.01,
         exclude=albens_common_para | sideroxylon_common_para |
          (albens_common_hwe01 & sideroxylon_common_hwe01) |
          is.na(albens_common_para) | is.na(sideroxylon_common_para))

test=myframe %>%
  filter(grepl("^Chr",snpID,fixed=F)) %>%
  mutate(chromo=factor(chromo),
         albens_common_para=factor(albens_common_para),
         albens_common_hwe01=factor(albens_common_hwe01),
         sideroxylon_common_para=factor(sideroxylon_common_para),
         sideroxylon_common_hwe01=factor(sideroxylon_common_hwe01),
         exclude=factor(exclude)
         ) %>%
  head(n=100000L) %>%
  # head() %>%
  tableplot(.,select=c(chromo,albens_common.hwe.Ho,albens_common_para,albens_common_hwe01,sideroxylon_common_para,sideroxylon_common_hwe01),
            sortCol = albens_common.hwe.Ho,nBins=100)

tableplot(test,select=c(chromo,albens_common.hwe.Ho,albens_common_para,albens_common_hwe01,sideroxylon_common_para,sideroxylon_common_hwe01),
          sortCol = albens_common.hwe.Ho,nBins=100)

tableplot(test)

head(test$columns$feral_hwe01)
