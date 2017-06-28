library(sp)
library(rgdal)

arca.dma <- readOGR('//deqhq1/tmdl/tmdl_wr/midcoast/Models/Sediment/SSN//lsn05_arca_sedstressor_03152017_dma.shp')

head(arca.dma@data)

library(dplyr)

dma_distinct <- arca.dma@data %>% group_by(rid_LSN05) %>% distinct(DMA_RP)

dma_distinct$rid_LSN05 <- as.factor(dma_distinct$rid_LSN05)

# Indian Creek DMA's
#        DMA_RP rid_LSN05
#         <chr>    <fctr>
# 1    USA-USFS     11495
# 2 ODF-PRIVATE     11495
# 3         ODA     11495
# 4 Lane County     11495'