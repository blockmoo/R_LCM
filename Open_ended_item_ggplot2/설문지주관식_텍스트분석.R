##############################################
library(wordcloud2)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(tidyverse)
# setwd("D:/myfriands10_google/삼성_설문지 분석/주관식/result")
setwd("C:/Users/hyerin/Google 드라이브(myfriands10@gmail.com)//삼성_설문지 분석/주관식/result")

dt1 <- read.csv("본인_19_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt2 <- read.csv("본인_20_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt3 <- read.csv("상사_19_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt4 <- read.csv("상사_20_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt5 <- read.csv("부서원_19_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt6 <- read.csv("부서원_20_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt7 <- read.csv("동료_19_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt8 <- read.csv("동료_20_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")

dt9 <- read.csv("본인_19_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt10 <- read.csv("본인_20_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt11 <- read.csv("상사_19_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt12 <- read.csv("상사_20_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt13 <- read.csv("부서원_19_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt14 <- read.csv("부서원_20_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt15 <- read.csv("동료_19_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt16 <- read.csv("동료_20_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")

dt17 <- read.csv("본인_19_기타.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt18 <- read.csv("상사_19_기타.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt19 <- read.csv("부서원_19_기타.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt20 <- read.csv("동료_19_기타.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")

dt21 <- read.csv("본인_21_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt22 <- read.csv("상사_21_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt23 <- read.csv("부서원_21_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt24 <- read.csv("동료_21_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")

dt25 <- read.csv("본인_21_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt26 <- read.csv("상사_21_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt27 <- read.csv("부서원_21_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt28 <- read.csv("동료_21_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")

dt29 <- read.csv("타인_19_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt30 <- read.csv("타인_20_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt31 <- read.csv("타인_21_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt32 <- read.csv("타인_19_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt33 <- read.csv("타인_20_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt34 <- read.csv("타인_21_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")

dt35 <- read.csv("임원본인_19_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt36 <- read.csv("임원본인_20_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt37 <- read.csv("임원본인_21_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt38 <- read.csv("임원본인_19_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt39 <- read.csv("임원본인_20_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt40 <- read.csv("임원본인_21_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")

dt41 <- read.csv("임원타인_19_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt42 <- read.csv("임원타인_20_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt43 <- read.csv("임원타인_21_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt44 <- read.csv("임원타인_19_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt45 <- read.csv("임원타인_20_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt46 <- read.csv("임원타인_21_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")

dt47 <- read.csv("비임원본인_19_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt48 <- read.csv("비임원본인_20_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt49 <- read.csv("비임원본인_21_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt50 <- read.csv("비임원본인_19_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt51 <- read.csv("비임원본인_20_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt52 <- read.csv("비임원본인_21_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")

dt53 <- read.csv("비임원타인_19_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt54 <- read.csv("비임원타인_20_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt55 <- read.csv("비임원타인_21_강점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt56 <- read.csv("비임원타인_19_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt57 <- read.csv("비임원타인_20_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")
dt58 <- read.csv("비임원타인_21_단점.csv", header = F,fileEncoding = "UTF-8", encoding = "euc-kr")


###

getwd()
#dir()

BOSS_DATA <- list()
boss_data_nm <- dir()[str_starts(dir(), pattern = "임원")] #디렉토리에서 임원으로 시작하는 데이터 불러오기
for(i in 1:length(boss_data_nm)) {
  data_sub <- read_csv(boss_data_nm[i], col_names = FALSE)
  BOSS_DATA[[i]] <- data_sub
}

NO_BOSS_DATA <- list()
no_boss_data_nm <- dir()[str_starts(dir(), pattern = "비임원")]
for(i in 1:length(no_boss_data_nm)) {
  data_sub <- read_csv(no_boss_data_nm[i], col_names = FALSE)
  NO_BOSS_DATA[[i]] <- data_sub
}

merged_data <- c(BOSS_DATA, NO_BOSS_DATA)

ex <- lapply(merged_data, FUN = function(x) mutate(x, total = sum(X2)) %>% 
               mutate(prob = round((X2/total)*100,2)))


names(ex) <- c(boss_data_nm, no_boss_data_nm)

names(ex)
# [1] "임원본인_19_강점.csv"   "임원본인_19_단점.csv"   "임원본인_20_강점.csv"  
# [4] "임원본인_20_단점.csv"   "임원본인_21_강점.csv"   "임원본인_21_단점.csv"  
# [7] "임원타인_19_강점.csv"   "임원타인_19_단점.csv"   "임원타인_20_강점.csv"  
# [10] "임원타인_20_단점.csv"   "임원타인_21_강점.csv"   "임원타인_21_단점.csv"  
# [13] "비임원본인_19_강점.csv" "비임원본인_19_단점.csv" "비임원본인_20_강점.csv"
# [16] "비임원본인_20_단점.csv" "비임원본인_21_강점.csv" "비임원본인_21_단점.csv"
# [19] "비임원타인_19_강점.csv" "비임원타인_19_단점.csv" "비임원타인_20_강점.csv"
# [22] "비임원타인_20_단점.csv" "비임원타인_21_강점.csv" "비임원타인_21_단점.csv"
####





##### ggplot
ggplot(dt1[1:10,] %>% mutate(color = ifelse(V1 == "경험", "blue", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 단어", "변화폭이 큰 단어"),values = c("black", "blue"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 본인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())

ggsave("본인_19_강점.png")
#-------------------------------------------
ggplot(dt2[1:10,]  %>% mutate(color = ifelse(V1 == "다양", "red",ifelse(V1 %in% c("경험","문제"),"blue", "black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점", "추가된 강점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 본인 top10 강점 ')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("본인_20_강점.png")

#-------------------------------------------
ggplot(dt21[1:10,] %>% mutate(color = ifelse(V1 == "추진", "red", ifelse(V1=="문제","blue","black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점" ,"추가된 강점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 본인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("본인_21_강점.png")
#-------------------------------------------

ggplot(dt3[1:10,] %>% mutate(color = ifelse(V1 == "지식", "blue", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점", "변화폭이 큰 강점"),values = c("black", "blue"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 상사 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("상사_19_강점.png")

#-------------------------------------------

ggplot(dt4[1:10,] %>% mutate(color = ifelse(V1 %in% c("전문성","책임감"), "red", ifelse(V1 =="지식","blue", "black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점", "추가된 강점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 상사 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("상사_20_강점.png")
#-------------------------------------------
ggplot(dt22[1:10,] %>% mutate(color = ifelse(V1 == "소통", "black", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점", "추가된 강점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 상사 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("상사_21_강점.png")
#-------------------------------------------
ggplot(dt5[1:10,] %>% mutate(color = ifelse(V1 == "능력", "blue", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점", "변화폭이 큰 강점"),values = c("black", "blue"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 부서원 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("부서원_19_강점.png")
#-------------------------------------------

ggplot(dt6[1:10,] %>% mutate(color = ifelse(V1 %in% c("지식","조직"), "red",ifelse(V1 == "능력","blue", "black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점", "추가된 강점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 부서원 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("부서원_20_강점.png")
#-------------------------------------------
ggplot(dt23[1:10,] %>% mutate(color = ifelse(V1 == "추진", "red", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점", "추가된 강점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 부서원 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("부서원_21_강점.png")
#-------------------------------------------

ggplot(dt7[1:10,] %>% mutate(color = ifelse(V1 == "지식", "blue", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점", "변화폭이 큰 강점"),values = c("black", "blue"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 동료 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("동료_19_강점.png")
#-------------------------------------------

ggplot(dt8[1:10,] %>% mutate(color = ifelse(V1 %in% c("협업","목표"), "red", ifelse(V1 == "지식","blue","black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점", "추가된 강점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 동료 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("동료_20_강점.png")
#-------------------------------------------
ggplot(dt24[1:10,] %>% mutate(color = ifelse(V1 %in% c("추진","역량"), "red", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점", "추가된 강점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 동료 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("동료_21_강점.png")
#-------------------------------------------

ggplot(dt9[1:10,] %>% mutate(color = ifelse(V1 %in% c("글로벌"), "purple", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점", "일년뒤 좋아진 약점"),values = c("black", "purple"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 본인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("본인_19_단점.png")
#-------------------------------------------

ggplot(dt10[1:10,] %>% mutate(color = ifelse(V1 %in% c("지식","학습"), "red",ifelse(V1 == "어학","purple", "black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점", "일년뒤 좋아진 약점","추가된 약점"),values = c("black","purple", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 본인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("본인_20_단점.png")
#-------------------------------------------
ggplot(dt25[1:10,] %>% mutate(color = ifelse(V1 %in% c("관리","개발"), "red", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점", "추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 본인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("본인_21_단점.png")
#-------------------------------------------

ggplot(dt11[1:10,] %>% mutate(color = ifelse(V1 %in% c("향상","추진"), "blue", ifelse(V1=="목표", "purple","black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점", "변화폭이 큰 약점","일년뒤 좋아진 약점"),values = c("black", "blue","purple"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 상사 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("상사_19_단점.png")
#-------------------------------------------

ggplot(dt12[1:10,] %>% mutate(color = ifelse(V1 == "개발", "red",ifelse(V1 %in% c("향상","추진"), "blue", "black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점", "변화폭이 큰 약점", "추가된 약점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 상사 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("상사_20_단점.png")
#-------------------------------------------
ggplot(dt26[1:10,] %>% mutate(color = ifelse(V1 =="가가", "blue", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점", "변화폭이 큰 약점"),values = c("black", "blue"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 상사 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("상사_21_단점.png")

#-------------------------------------------

ggplot(dt13[1:10,] %>% mutate(color = ifelse(V1 %in% c("기술"),"blue", ifelse(V1=="관심","purple", "black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점", "일년뒤 좋아진 약점"),values = c("black","blue", "purple"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 부서원 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("부서원_19_단점.png")
#-------------------------------------------

ggplot(dt14[1:10,] %>% mutate(color = ifelse(V1 %in% c("역량","방향"), "red", ifelse(V1 %in% c("기술"),"blue",ifelse(V1== "지시","purple","black")))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점", "일년뒤 좋아진 약점", "추가된 약점"),values = c("black","blue","purple", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 부서원 top10 약점 ')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("부서원_20_단점.png")
#-------------------------------------------
ggplot(dt27[1:10,] %>% mutate(color = ifelse(V1 %in% c("건강","그룹","지금","관심"), "red","black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black","red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 부서원 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("부서원_21_단점.png")
#-------------------------------------------

ggplot(dt15[1:10,] %>% mutate(color = ifelse(V1 %in% c("협업","능력"),"blue","black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점"),values = c("black", "blue"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 동료 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("동료_19_단점.png")
#-------------------------------------------

ggplot(dt16[1:10,] %>% mutate(color = ifelse(V1 == "리더", "red", ifelse(V1 %in% c("협업","능력","관리","의견"),"blue","black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점", "추가된 약점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 동료 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("동료_20_단점.png")

#-------------------------------------------
ggplot(dt28[1:10,] %>% mutate(color = ifelse(V1 %in% c("후배","건강"), "red", ifelse(V1 %in% c("관리","의견"),"blue","black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점", "추가된 약점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 동료 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("동료_21_단점.png")


#-------------------------------------------

ggplot(dt17[1:10,] %>% mutate(color = ifelse(V1 == "소통", "red", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 단점", "추가된 단점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 본인 기타 상위 10개')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("본인_19_기타.png")
#-------------------------------------------

ggplot(dt18[1:10,] %>% mutate(color = ifelse(V1 == "소통", "red", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 단점", "추가된 단점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 상사 기타 상위 10개')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("상사_19_기타.png")
#-------------------------------------------

ggplot(dt19[1:10,] %>% mutate(color = ifelse(V1 == "소통", "red", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 단점", "추가된 단점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 부서원 기타 상위 10개')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("부서원_19_기타.png")
#-------------------------------------------

ggplot(dt20[1:10,] %>% mutate(color = ifelse(V1 == "소통", "red", "black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 단점", "추가된 단점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 동료 기타 상위 10개')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("동료_19_기타.png")

#--------------------------------------------
#--------------------------------------------
ggplot(dt29[1:10,] %>% mutate(color = ifelse(V1 %in% c("추진","의견"), "blue","black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점"),values = c("black","blue"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 타인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("타인_19_강점.png")


#-------------------------------------------
ggplot(dt30[1:10,] %>% mutate(color = ifelse(V1 %in% c("지식"), "red", ifelse(V1 %in% c("추진","의견", "소통"),"blue","black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점", "추가된 약점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 타인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("타인_20_강점.png")


#-------------------------------------------
ggplot(dt31[1:10,] %>% mutate(color = ifelse(V1 %in% c("조직"), "red", ifelse(V1 %in% c("추진","소통"),"blue","black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점", "추가된 약점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 타인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("타인_21_강점.png")


#-------------------------------------------
ggplot(dt32[1:10,] %>% mutate(color = ifelse(V1 %in% c("조직"), "blue", ifelse(V1 %in% c("관심"),"purple","black"))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점", "일년뒤 좋아진 약점"),values = c("black","blue","purple"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 타인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("타인_19_단점.png")


#-------------------------------------------
ggplot(dt33[1:10,] %>% mutate(color = ifelse(V1 %in% c("방향"), "red", ifelse(V1 %in% c("조직"),"blue",ifelse(V1 %in% c("능력","적극","추진"),"purple","black")))), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점","일년뒤 좋아진 약점", "추가된 약점"),values = c("black","blue", "purple","red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 타인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("타인_20_단점.png")


#-------------------------------------------
ggplot(dt34[1:10,] %>% mutate(color = ifelse(V1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(V1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 타인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("타인_21_단점.png")


#-------------------------------------------

#########################################################
###### 임원과 비임원 비교

### 임원본인 강점
ggplot(ex[[1]][1:10,] %>% mutate(color = ifelse(X1 %in% c("추진","경험"), "blue","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점"),values = c("black", "blue"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 임원본인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원본인_19_강점.png")
#-------------------------------------------
ggplot(ex[[3]][1:10,] %>% mutate(color = ifelse(X1 %in% c("역량","결정"), "red",ifelse(X1 %in% c("경험","추진"),"blue","black"))), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점","추가된 강점"),values = c("black","blue", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 임원본인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원본인_20_강점.png")
#-------------------------------------------
ggplot(ex[[5]][1:10,] %>% mutate(color = ifelse(X1 %in% c("능력","의견"), "red",ifelse(X1 %in% c("추진","해결"),"blue","black"))), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점","추가된 강점"),values = c("black", "blue","red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 임원본인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원본인_21_강점.png")

### 임원본인 단점
#-------------------------------------------
ggplot(ex[[2]][1:10,] %>% mutate(color = ifelse(X1 %in% c("지속","추진","성과","이해"), "purple",ifelse(X1 %in% c("소통"),"blue","black"))), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점","일년뒤 좋아진 약점"),values = c("black", "blue","purple"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 임원본인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원본인_19_단점.png")
#-------------------------------------------
ggplot(ex[[4]][1:10,] %>% mutate(color = ifelse(X1 %in% c("지식","조직","분야","관리"), "red",ifelse(X1 %in% c("소통"),"blue",ifelse(X1 %in% c("사업","글로벌"),"purple","black")))), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점",  "일년뒤 좋아진 약점", "추가된 약점"),values = c("black","blue","purple", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 임원본인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원본인_20_단점.png")
#-------------------------------------------
ggplot(ex[[6]][1:10,] %>% mutate(color = ifelse(X1 %in% c("추진","결정","이해","경영"), "red",ifelse(X1 %in% c("소통"),"blue",
"black"))), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","변화폭이 큰 약점","추가된 약점"),values = c("black", 'blue',"red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 임원본인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원본인_21_단점.png")
#-------------------------------------------

### 임원타인 강점
ggplot(ex[[7]][1:10,] %>% mutate(color = ifelse(X1 %in% c("해결","문제"), "blue","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점"),values = c("black", "blue"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 임원타인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원타인_19_강점.png")
#-------------------------------------------
ggplot(ex[[9]][1:10,] %>% mutate(color = ifelse(X1 %in% c("방향"), "red",ifelse(X1 %in% c("해결","문제", "경험","소통"),"blue","black"))), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점","추가된 강점"),values = c("black", "blue","red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 임원타인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원타인_20_강점.png")
#-------------------------------------------
ggplot(ex[[11]][1:10,] %>% mutate(color = ifelse(X1 %in% c("제시"), "red",ifelse(X1 %in% c("소통","방향","경험"),"blue","black"))), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 강점","변화폭이 큰 강점","추가된 강점"),values = c("black","blue" ,"red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 임원타인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원타인_21_강점.png")
#-------------------------------------------

### 임원타인 단점
ggplot(ex[[8]][1:10,] %>% mutate(color = ifelse(X1 %in% c("이해","능력","적극"), "purple","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","일년뒤 좋아진 약점"),values = c("black", "purple"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 임원타인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원타인_19_단점.png")
#-------------------------------------------
ggplot(ex[[10]][1:10,] %>% mutate(color = ifelse(X1 %in% c("관리","결정","방향"), "red",ifelse(X1 %in% c("역량","사업"),"purple","black"))), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","일년뒤 좋아진 약점","추가된 약점"),values = c("black", "purple","red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 임원타인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원타인_20_단점.png")
#-------------------------------------------
ggplot(ex[[12]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","기회","직원"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 임원타인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("임원타인_21_단점.png")
#-------------------------------------------

### 비임원본인 강점
ggplot(ex[[13]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 비임원본인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원본인_19_강점.png")
#-------------------------------------------
ggplot(ex[[15]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 비임원본인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원본인_20_강점.png")
#-------------------------------------------
ggplot(ex[[17]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 비임원본인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원본인_21_강점.png")
#-------------------------------------------

### 비임원본인 단점
ggplot(ex[[14]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 비임원본인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원본인_19_단점.png")
#-------------------------------------------
ggplot(ex[[16]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 비임원본인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원본인_20_단점.png")
#-------------------------------------------
ggplot(ex[[18]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 비임원본인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원본인_21_단점.png")
#-------------------------------------------

### 비임원타인 강점
ggplot(ex[[19]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 비임원타인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원타인_19_강점.png")
#-------------------------------------------
ggplot(ex[[21]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 비임원타인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원타인_20_강점.png")
#-------------------------------------------
ggplot(ex[[23]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 비임원타인 top10 강점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원타인_21_강점.png")
#-------------------------------------------
ggplot(ex[[20]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2019 비임원타인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원타인_19_단점.png")
#-------------------------------------------
ggplot(ex[[22]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2020 비임원타인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원타인_20_단점.png")
#-------------------------------------------
ggplot(ex[[24]][1:10,] %>% mutate(color = ifelse(X1 %in% c("건강","관심","그룹"), "red","black")), aes(prob, reorder(X1,prob), fill = color)) +
  geom_col() + 
  geom_text(aes(x=prob ,label=paste0(prob,"%")),colour='white', 
            position = position_stack(vjust=0.8),size=6) +
  theme(axis.text.y=element_text(size=20)) +
  scale_fill_manual(labels=c("기존 약점","추가된 약점"),values = c("black", "red"))+
  labs(x = '단어 비율(단어 빈도수/전체 단어 빈도수)',
       y = '',
       title = '2021 비임원타인 top10 약점')+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 15),
        title = element_text(size=20),
        panel.background = element_blank())
ggsave("비임원타인_21_단점.png")
#-------------------------------------------




























