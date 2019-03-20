

# devtools::install_github("matthewkling/chelsaDL")
# devtools::install_github("matthewkling/climatica")

library(chelsaDL)
library(climatica)
library(tidyverse)
library(raster)

select <- dplyr::select



############ download monthly climate variables from CHELSA ##########################

# study area bounding box
# (for cropping climate data, to limit file sizes and processing times)
fia <- read_csv("data/fia.csv") %>%
      filter(lon < -95, lat < 50, lat > 25) %>%
      select(gs, lat, lon)
coordinates(fia) <- c("lon", "lat")
bbox <- extent(fia)

# define variables and download
ch <- ch_queries(variables=c("tmin", "tmax", "temp", "prec"),
                 models=c("HadGEM2-ES", "CNRM-CM5", "CanESM2", "MIROC5"),
                 scenarios=c("rcp45", "rcp85"),
                 timeframes=c("2061-2080"),
                 months=1:12) %>%
      ch_dl(dest="f:/chelsa/toporefugia/future/raw_monthly",
            crop=bbox)


# one model is missing mean temperature. 
# derive it from tmin and tmax.
missing <- ch %>% filter(status == "incomplete", variable=="temp")
if(nrow(missing) > 0){
      for(i in 1:nrow(missing)){
            m <- missing[i,]
            n <- ch %>% filter(model==m$model, scenario==m$scenario, 
                               month==m$month, timeframe==m$timeframe,
                               variable %in% c("tmin", "tmax")) %>%
                  pull(path)
            
            # build outfile name, with correct model run number
            out <- str_split(n[1], "i1p1")[[1]][1] %>% 
                  substr(nchar(.)-2, nchar(.)) %>%
                  str_extract("\\d") %>%
                  sub("\\*", ., m$path)
            
            # derive mean temperature
            n %>% stack() %>%
                  mean() %>%
                  writeRaster(out)
            
            # update the metadata list
            chi <- ch$model==m$model & ch$scenario==m$scenario & 
                  ch$month==m$month & ch$timeframe==m$timeframe &
                  ch$variable=="temp"
            ch$path[chi] <- out
            ch$status[chi] <- "derived"
      }
}

# historic data -- already downloaded, just crop and 
for(file in list.files("F:/chelsa/monthly48", full.names=T)){
      out <- sub("\\.tif", ".grd", basename(file)) %>%
            paste0("f:/chelsa/toporefugia/historic/raw_monthly/", .)
      file %>% raster() %>% crop(bbox) %>% writeRaster(out)
}


############# derive new climate variables ################

derive <- function(md, outdir){
      # md <- batches[[2]]
      outfile <- paste0(outdir, md$set[1], ".grd")
      if(file.exists(outfile)) return("skipped")
      message(paste("processing", md$set[1]))
      
      # load layers, in correct order for water balance function
      md <- md %>% 
            mutate(var_ord=case_when(variable=="prec" ~ 1,
                                     variable=="temp" ~ 2,
                                     variable=="tmax" ~ 3,
                                     variable=="tmin" ~ 4),
                   month=str_pad(month, 2, "left", 0)) %>%
            arrange(var_ord, month)
      r <- stack(md$path)
      
      # water balance variables
      water <- hydro(r, temp_scalar=0.1, ncores=6, already_latlong=T)
      names(water) <- tolower(names(water))
      
      # temprature variables
      temps <- stack()
      temps$mat <- md %>% filter(variable=="temp") %>%
            pull(path) %>%  stack() %>% mean() %>% "*"(0.1)
      temps$tmax <- md %>% filter(variable=="tmax") %>%
            pull(path) %>%  stack() %>% mean() %>% "*"(0.1)
      temps$tmin <- md %>% filter(variable=="tmin") %>%
            pull(path) %>%  stack() %>% mean() %>% "*"(0.1)
      temps$djf <- md %>% filter(variable=="tmin", 
                                 month %in% c("12", "01", "02")) %>%
            pull(path) %>%  stack() %>% mean() %>% "*"(0.1)
      temps$jja <- md %>% filter(variable=="tmax", 
                                 month %in% c("06", "06", "08")) %>%
            pull(path) %>%  stack() %>% mean() %>% "*"(0.1)
      temps$tminmin <- md %>% filter(variable=="tmin") %>%
            pull(path) %>%  stack() %>% min() %>% "*"(0.1)
      temps$tmaxmax <- md %>% filter(variable=="tmax") %>%
            pull(path) %>%  stack() %>% max() %>% "*"(0.1)
      
      # export
      writeRaster(stack(water, temps), outfile, overwrite=T)
      return(outfile)
}


ch %>% 
      mutate(set=paste(timeframe, scenario, model, sep="_")) %>%
      arrange(set) %>%
      split(.$set) %>%
      lapply(derive, outdir="f:/chelsa/toporefugia/future/derived/")

parseMetadata("f:/chelsa/toporefugia/historic/raw_monthly/", pattern="\\.grd",
              keys=list(var=c("tmax10", "tmin10", "temp10", "prec"), mo=1:12)) %>%
      select(-month) %>%
      rename(variable=var, month=mo) %>%
      mutate(variable=sub(10, "", variable),
             set="historic") %>%
      derive(outdir="f:/chelsa/toporefugia/historic/derived/")


