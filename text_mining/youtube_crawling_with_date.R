setwd("C:/Users/hyerin/Google 드라이브(myfriands10@gmail.com)/Mixture_lab/youtube_산출물/R결과물2")
setwd("D:/myfriands10_google/Mixture_lab/youtube_산출물/R결과물2")
setwd("/Users/hagsaeng/Downloads/selenium")
setwd("C:/Users/hyerin/Google 드라이브(myfriands10@gmail.com)/Mixture_lab/youtube_산출물/댓글재수집")


kid <- read.csv('date_k_kids_disaster.csv')
covid <- read.csv("date_korea_covid19_kids.csv")
u.kid <- read.csv("date_usa_kids_disaster.csv")
u.covid <- read.csv("date_usa_covid19_kids.csv")

dim(kid) #11929
dim(covid) #13403
dim(u.kid) #17610
dim(u.covid) #47726


### **Keywords **

키워드 | 동영상 갯수 | 댓글수
-|-|-
코로나19 + 유아 |549 개  |17588 개
유아 + 재난안전 |518 개  |7819 개
COVID-19 + children |570 개  |38009 개
children + disaster safety |458 개  |79419 개

* 날짜 추가 후 (04/08)
키워드 | 동영상 갯수 | 댓글수
-|-|-
코로나19 + 유아 |549 개  |10978 개
유아 + 재난안전 |518 개  |5874 개
COVID-19 + children |570 개  |12624 개
children + disaster safety |458 개  |55288 개

* 날짜 추가 후 (05/06)
키워드 | 동영상 갯수 | 댓글수
-|-|-
코로나19 + 유아 |549 개  |15745 개
유아 + 재난안전 |518 개  |2436 개
COVID-19 + children |570 개  |18135 개
children + disaster safety |458 개  |39084 개




# 04/22 교수님 컴으로 재크롤링한 결과
dim(new.u.disaster) #55088     4
dim(new.u.covid) #10330     4
dim(new.k.covid) #8100    4
dim(new.k.disaster)  #5700    4


library(dplyr)
library(tidyverse)
library(stringr)

#### extract links that can be used
arr_day <-function(x){
  x$day <- str_replace_all(x$url_day,"실시간 스트리밍 시작일: ","")
  x$day <- str_replace_all(x$url_day,"최초 공개:","")
  x$day <- str_trim(x$url_day)
  x <- x[which(!x$url_day %in% c("비공개 동영상","동영상을 재생할 수 없음")),]
  
  x <- x %>%
    separate(url_day, into = c("year", "month", "day")) %>%
    unite(col = "date", year, month, day, sep = "-") %>%
    mutate(date = as.Date(date)) %>%
    arrange(date)
  return(x)
}

covid.child <- arr_day(u.covid)
k.dis <- arr_day(kid)
k.covid.child <- arr_day(covid)
dis <- arr_day(u.kid)
dim(kid) #11929
dim(covid) #13403
dim(u.kid) #17610
dim(u.covid) #47726

View(dis)

####################################################################################
###### extract ~2020.12.29 comments
setwd("C:/Users/hyerin/Google 드라이브(myfriands10@gmail.com)/Mixture_lab/youtube_산출물/댓글재수집")
setwd("D:/myfriands10_google/Mixture_lab/youtube_산출물/댓글재수집")

library(stringr)

k.disaster <- read.csv("youtube_k_kids_disaster.csv")
k.covid <- read.csv("youtube_k_covid19_kids.csv")
u.kid.covid <- read.csv("youtube_usa_covid19_kid.csv")
u.kid.disaster <- read.csv("youtube_usa_kids_disaster.csv")
k <- read.csv("youtube_korea.csv")

View(k.disaster)
# 1.
table(k.disaster$report_date)
stop <- c("일 전", "주 전", "시간 전", "1개월", "2개월", "3개월")
new.k.disaster <- k.disaster
new.k.disaster <- new.k.disaster[-grep("시간 전", new.k.disaster$report_date),]
new.k.disaster <- new.k.disaster[-grep("일 전", new.k.disaster$report_date),]
new.k.disaster <- new.k.disaster[-grep("주 전", new.k.disaster$report_date),]
new.k.disaster <- new.k.disaster[-grep("1개월", new.k.disaster$report_date),]
new.k.disaster <- new.k.disaster[-grep("2개월", new.k.disaster$report_date),]
new.k.disaster <- new.k.disaster[-grep("3개월", new.k.disaster$report_date),]
new.k.disaster$report_date <- gsub("\\(수정됨)","",new.k.disaster$report_date)
table(new.k.disaster$report_date)
dim(new.k.disaster)

# 2.
table(k.covid$report_date)
dim(k.covid)
stop <- c("일 전", "주 전", "시간 전", "1개월", "2개월", "3개월", "분 전")
data <- k.covid
data <- data[-grep("시간 전", data$report_date),]
data <- data[-grep("일 전", data$report_date),]
data <- data[-grep("주 전", data$report_date),]
data <- data[-grep("1개월", data$report_date),]
data <- data[-grep("2개월", data$report_date),]
data <- data[-grep("3개월", data$report_date),]
data <- data[-grep("분 전", data$report_date),]
data$report_date <- gsub("\\(수정됨)","",data$report_date)
new.k.covid <- data

table(new.k.covid$report_date)
dim(new.k.covid)


# 3.
table(u.kid.covid$report_date)
dim(u.kid.covid)
stop <- c("일 전", "주 전", "시간 전", "1개월", "2개월", "3개월", "분 전")
data <- u.kid.covid
data <- data[-grep("시간 전", data$report_date),]
data <- data[-grep("일 전", data$report_date),]
data <- data[-grep("주 전", data$report_date),]
data <- data[-grep("1개월", data$report_date),]
data <- data[-grep("2개월", data$report_date),]
data <- data[-grep("3개월", data$report_date),]
data <- data[-grep("분 전", data$report_date),]
data$report_date <- gsub("\\(수정됨)","",data$report_date)
new.u.covid <- data

table(new.u.covid$report_date)
dim(new.u.covid)


# 4.
table(u.kid.disaster$report_date)
dim(u.kid.disaster)
stop <- c("일 전", "주 전", "시간 전", "1개월", "2개월", "3개월", "분 전")
data <- u.kid.disaster
data <- data[-grep("시간 전", data$report_date),]
data <- data[-grep("일 전", data$report_date),]
data <- data[-grep("주 전", data$report_date),]
data <- data[-grep("1개월", data$report_date),]
data <- data[-grep("2개월", data$report_date),]
data <- data[-grep("3개월", data$report_date),]
data <- data[-grep("분 전", data$report_date),]
data$report_date <- gsub("\\(수정됨)","",data$report_date)
new.u.disaster <- data
table(new.u.disaster$report_date)
dim(new.u.disaster)


dim(new.k.disaster) #5874
dim(new.k.covid) #10978
dim(new.u.covid) #12624
dim(new.u.disaster) #55288
names(new.k.disaster)
#
new.k.dis <- new.k.disaster %>% arrange(report_date)
new.k.dis <- new.k.dis %>%   group_by(report_date) %>% 
  summarise(n=n()) %>% as.data.frame()
new.k.dis
#
new.k.co <- new.k.covid %>% arrange(report_date)
new.k.co <- new.k.co %>%   group_by(report_date) %>% 
  summarise(n=n()) %>% as.data.frame()
new.k.co
#
new.u.co <- new.u.covid %>% arrange(report_date)
new.u.co <- new.u.co %>%   group_by(report_date) %>% 
  summarise(n=n()) %>% as.data.frame()
new.u.co
#
new.u.dis <- new.u.disaster %>% arrange(report_date)
new.u.dis <- new.u.dis %>%   group_by(report_date) %>% 
  summarise(n=n()) %>% as.data.frame()
new.u.dis


save(new.k.co, new.k.dis, new.u.co, new.u.dis,file = "./new_comments_num.r")



#################################
###### 05/06까지 수집한 결과 확인
setwd("D:/myfriands10_google/Mixture_lab/youtube_산출물/R결과물2/날짜추가댓글_수정본")

k.disaster <- read.csv("0503youtube_k_kids_disaster.csv")
k.covid <- read.csv("0504youtube_korea_covid19_kids.csv")
u.kid.covid <- read.csv("0429youtube_usa_covid19_kid.csv")
u.kid.disaster <- read.csv("0504youtube_usa_kids_disaster.csv")


new.k.dis <- read.csv("new.k.disaster.csv")
new.k.cov <- read.csv("new.k.covid.csv")
new.u.dis <- read.csv("new.u.kid.covid.csv")
new.u.cov <- read.csv("new.u.kid.disaster.csv")


#############################
###### 05/06 3~20개의 topic num일때의 LDA결과를 list로 저장  
load("D:/myfriands10_google/Mixture_lab/youtube_산출물/R결과물/lda_list4.r")

# make colnames function : ex) topic2(5.8%)
find.pct <- function(x,y){
  top.pct <- c()
  for(i in 1:ncol(x)){
    top.pct <- c(top.pct,paste0("주제 ",i,"(",y$topic.pct[i],")"))
  }
  return(top.pct)
}

# save results : topic num(3~20)
df <- lapply(k.dis.list, FUN = function(x) as.data.frame(x$new_LDA_result) )
name.list <- mapply(find.pct, x=df,y=k.dis.list)
t.k.dis <- mapply(FUN = function(x,y) {names(x) <- y; return(x)}, x=df,y=name.list)
# t.k.dis <- mapply(FUN = function(x,y) {names(x) <- y$topic.pct; return(x)}, x=df,y=k.dis.list)

df <- lapply(k.covid.list, FUN = function(x) as.data.frame(x$new_LDA_result) )
name.list <- mapply(find.pct, x=df,y=k.covid.list)
t.k.cov <- mapply(FUN = function(x,y) {names(x) <- y; return(x)}, x=df,y=name.list)
# t.k.cov <- mapply(FUN = function(x,y) {names(x) <- y$topic.pct; return(x)}, x=df,y=k.covid.list)

df <- lapply(u.dis.list, FUN = function(x) as.data.frame(x$new_LDA_result) )
name.list <- mapply(find.pct, x=df,y=u.dis.list)
t.u.dis <- mapply(FUN = function(x,y) {names(x) <- y; return(x)}, x=df,y=name.list)
# t.u.dis <- mapply(FUN = function(x,y) {names(x) <- y$topic.pct; return(x)}, x=df,y=u.dis.list)

df <- lapply(u.covid.list, FUN = function(x) as.data.frame(x$new_LDA_result) )
name.list <- mapply(find.pct, x=df,y=u.covid.list)
t.u.cov <- mapply(FUN = function(x,y) {names(x) <- y; return(x)}, x=df,y=name.list)


lda_list4 <- list(t.k.cov, t.k.dis, t.u.cov, t.u.dis)

save(t.k.dis, t.k.cov, t.u.dis, t.u.cov,lda_list4, file = "./R결과물/table_lda_list4.r")
load("./R결과물/table_lda_list4.r")

# -> 결과는 shiny를 이용해서 볼 수 있도록 정리.
# ./R결과물/LDA_topic_num의 app.R
# https://youtube-hyerin.shinyapps.io/LDA_topic_num_hyerin/


