library(officer)
library(flextable)
library(dplyr)
library(purrr)
library(readxl)

find_slide <- function(txt) {
  tmp <- map(seq_along(pres), ~ {
    slide_summary(on_slide(pres, .x)) %>%
    mutate(slide_no = .x)
  }) 
  tmp2 <- tmp %>% map(~ {
    filter(.x, type=="title") %>%
    select(text, slide_no) %>% 
    mutate(text=unlist(text))
  }) %>% 
    bind_rows()
  tmp2$slide_no[which(tmp2$text==txt)]
}

path <- "D:\\Text\\Teaching\\NumericalCourse\\BES Workshop\\BES Workshop.pptx"

pres <- read_pptx(path)

layout_summary(pres)

pres <- add_slide(pres, layout="Title")

n <- find_slide("Tidy data")

pres <- on_slide(pres, 11)

hock <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet="HOCKHAM")

hock <- hock[1:5, c(1, 3, 4:6)]
ft <- flextable(hock)
ft <- autofit(ft)

loc1 <- ph_location(left=1, top=2, width=2, height=5, bg="lightgrey")

pres <- ph_with(pres, hock, location = loc1)


print(pres, "Test.pptx")
