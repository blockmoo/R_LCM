########################## packages

library(N2H4)
library(rvest)
library(rJava)
library(KoNLP)
library(data.table)
library(reshape)
library(dplyr)
library(wordcloud2)
library(RColorBrewer)
library(lda)
library(ldatuning)
library(topicmodels) #topicmodels, lda 
library(tm) # text mining, corpus, tdm, dim 등 제공. corpus형태로 저장. 
library(qdap) # 원 텍스트데이터를 dataframe으로 저장
library(stringr)
# library(wordcloud)
# library(SnowballC)  # for text stemming(어근추출)
library(htmlwidgets)
# install.packages("webshot")
library(webshot)
# webshot::install_phantomjs()
library(slam)
useNIADic() #1213109단어 built



########################## data preprocessing
# setwd("C:/Users/hyerin/Google 드라이브(myfriands10@gmail.com)/Mixture_lab/youtube_산출물")
setwd("D:/myfriands10_google/Mixture_lab/youtube_산출물")
child.dis <- read.csv("children_disaster_safety.csv")
covid.child <- read.csv("COVID_19_children.csv")
k.child.dis <- read.csv("유아_재난안전.csv")
k.covid.child <- read.csv("코로나19_유아.csv")

dim(child.dis)  #79419
dim(covid.child)  #38009
dim(k.covid.child) #17588
dim(k.child.dis) #7819



##################### 1. 유아+재난안전 
comments = readLines("유아_재난안전.csv") #readLines : 한줄씩 읽어 문자열 벡터로 저장. 
comments <- str_replace_all(comments, "\\<[^>]+\\>","") #이모티콘 제거
comments <- str_replace_all(comments, "\\W"," ") #특수문자 제거
head(comments)



#명사 추출 : 오류나는 문장들은 csv에서 수정 
nouns <- extractNoun(comments) 
length(nouns) #8606
head(nouns,4)


# 1차 스크리닝 (불필요한 단어들 삭제)  
for(i in 1:length(nouns)){
  nouns[[i]] <- gsub('[~!@#$%&*()_+=?<>]','',nouns[[i]])
  nouns[[i]] <- gsub("\\[","",nouns[[i]])
  nouns[[i]] <- gsub('[ㄱ-ㅎ]','',nouns[[i]])
  nouns[[i]] <- gsub('(ㅜ|ㅠ)','',nouns[[i]])
  nouns[[i]] <- gsub("\\^","",nouns[[i]])
  nouns[[i]] <- gsub("\\d","",nouns[[i]])
  nouns[[i]] <- gsub("[a-z]","",nouns[[i]])
  nouns[[i]] <- gsub("[A-Z]","",nouns[[i]])
  nouns[[i]] <- gsub("애","아이",nouns[[i]]) # 애 -> 아이
  nouns[[i]] <- gsub("아이들","아이",nouns[[i]]) # 아이들 -> 아이
  nouns[[i]] <- gsub("아가","아기",nouns[[i]]) # 아가-> 아기
  
}
head(nouns)
# View(nouns)

########################## word

stopwords = c("들이", "하게" ,"드가", "해서", "이거", "하시","하면", "도티",
              "하지", "저거", "하기", "해요", "해주", "하세")

### wordcloud
word_cloud <- function(nouns, stopwords){
  
  wordcount <- table(unlist(nouns))
  wordcount <- sort(wordcount, decreasing = TRUE)
  
  df_word <- as.data.frame(wordcount, stringsAsFactors = F)
  df_word$Var1=as.character(df_word$Var1)
  
  two <- filter(df_word, nchar(Var1) >= 2)
  one <- filter(df_word, nchar(Var1) == 1)
  one$Var1 <- gsub("\\d+","",one$Var1)
  one <- filter(one, one$Var1 != "") #NA 지우기
  # head(one,30)
  
  # 2차 스크리닝 
  stopword <- stopwords
  two_clean <- two %>% filter(!Var1 %in% stopword)
  # head(two_clean,30)
  
  # 한글자에서 의미있는 것 + 두글자 이상 단어들 
  # two_clean_kid <- rbind(two_clean, one[which(one$Var1 == "애"),]) %>% arrange(desc(Freq))
  # rownames(two_clean_kid) <- c(1:nrow(two_clean_kid))
  # head(two_clean_kid,30)
  
  #워드클라우드
  word.plot2 <- wordcloud2(data = two_clean)
  
  return(list(nouns=nouns, one=one, two = two,
              two_clean=two_clean, 
              wordcloud = word.plot2))
}

word <- word_cloud(nouns, stopwords)
head(word$one, 30)
head(word$two, 30)
head(word$two_clean, 30)
head(word$two_clean_kid, 30)
word.plot <- word$wordcloud


#html형식으로 워드클라우드 파일 저장하기
saveWidget(word.plot, "./result/tmp1.html", selfcontained = F)

#tmp.html 파일을 png로 변환하기
webshot("./result/tmp1.html", "./result/wordcloud1.png", delay =20, vwidth = 1200, vheight=900)


save(comments, nouns, stopwords, word, file = "./R결과물/k.disaster2.r")
load("./R결과물/k.disaster2.r")





##################### 2. 코로나19 + 유아  
comments = readLines("코로나19_유아.csv") #readLines : 한줄씩 읽어 문자열 벡터로 저장. 
comments <- str_replace_all(comments, "\\<[^>]+\\>","") #이모티콘 제거
comments <- str_replace_all(comments, "\\W"," ") #특수문자 제거
head(comments)



#명사 추출 : 오류나는 문장들은 csv에서 수정 
nouns <- extractNoun(comments) 
length(nouns) # 18252
head(nouns,4)


# 1차 스크리닝 (불필요한 단어들 삭제)  
for(i in 1:length(nouns)){
  nouns[[i]] <- gsub('[~!@#$%&*()_+=?<>]','',nouns[[i]])
  nouns[[i]] <- gsub("\\[","",nouns[[i]])
  nouns[[i]] <- gsub('[ㄱ-ㅎ]','',nouns[[i]])
  nouns[[i]] <- gsub('(ㅜ|ㅠ)','',nouns[[i]])
  nouns[[i]] <- gsub("\\^","",nouns[[i]])
  nouns[[i]] <- gsub("\\d","",nouns[[i]])
  nouns[[i]] <- gsub("[a-z]","",nouns[[i]])
  nouns[[i]] <- gsub("[A-Z]","",nouns[[i]])
  nouns[[i]] <- gsub("애","아이",nouns[[i]]) # 애 -> 아이
  nouns[[i]] <- gsub("아이들","아이",nouns[[i]]) # 아이들 -> 아이
  nouns[[i]] <- gsub("아가","아기",nouns[[i]]) # 아가-> 아기
}
head(nouns)
# View(nouns)

########################## word
#도티 추가함
stopwords = c("저희","들이", "하게" ,"드가", "해서", "이거", 
              "하시","하면", "들이", "이거", "때문", "로라",
              "에에에에에", "에에에ㅔ에에에에에", "도티",
              "하지", "저거", "하기", "해요", "해주")  #도티 추가함

### wordcloud
word_cloud <- function(nouns, stopwords){
  
  wordcount <- table(unlist(nouns))
  wordcount <- sort(wordcount, decreasing = TRUE)
  
  df_word <- as.data.frame(wordcount, stringsAsFactors = F)
  df_word$Var1=as.character(df_word$Var1)
  
  two <- filter(df_word, nchar(Var1) >= 2)
  one <- filter(df_word, nchar(Var1) == 1)
  one$Var1 <- gsub("\\d+","",one$Var1)
  one <- filter(one, one$Var1 != "") #NA 지우기
  # head(one,30)
  
  # 2차 스크리닝 
  stopword <- stopwords
  two_clean <- two %>% filter(!Var1 %in% stopword)
  # head(two_clean,30)
  
  # # 한글자에서 의미있는 것 + 두글자 이상 단어들 
  # two_clean_kid <- rbind(two_clean, one[which(one$Var1 == "애"),]) %>% arrange(desc(Freq))
  # rownames(two_clean_kid) <- c(1:nrow(two_clean_kid))
  # head(two_clean_kid,30)
  
  #워드클라우드
  word.plot2 <- wordcloud2(data = two_clean)
  
  return(list(nouns=nouns, one=one, two = two,
              two_clean=two_clean,
              wordcloud = word.plot2))
}


word <- word_cloud(nouns, stopwords)
head(word$one, 30)
head(word$two, 30)
head(word$two_clean, 30)
head(word$two_clean_kid, 30)
word$wordcloud
word.plot <- word$wordcloud


#html형식으로 워드클라우드 파일 저장하기
saveWidget(word.plot, "./result/tmp2.html", selfcontained = F)

#tmp.html 파일을 png로 변환하기
webshot("./result/tmp2.html", "./result/wordcloud2.png", delay =20, vwidth = 1200, vheight=900)


save(comments, nouns, stopwords, word, file = "./R결과물/k.covid2.r")
load("./R결과물/k.covid2.r")



########################### 3. COVID_19 + children

## 댓글 쪼개기 
covid <- read.csv("COVID_19_children.csv")
dim(covid) #38009 1 
covid1 <- covid[1:15000,]
write.csv(covid1, "COVID1.csv")
covid2 <- covid[15001:30000,]
write.csv(covid2, "COVID2.csv")
covid3 <- covid[30001:38009,]
write.csv(covid2, "COVID3.csv")


clean.cps <- function(comment){
  comments <- enc2utf8(comments)  # UTF-8로 변환
  cps <- Corpus(VectorSource(comments))
  cps.tm <- tm_map(cps, removePunctuation) # 구두점 제거
  cps.tm <- tm_map(cps.tm, removeNumbers) # 숫자제거
  cps.tm <- tm_map(cps.tm, tolower)       #알파벳 모두 소문자로 변경 
  cps.tm <- tm_map(cps.tm, stripWhitespace) # 공백 제거 
  cps.tm <- tm_map(cps.tm,removeWords, c(stopwords("english"),stopwords("SMART"),"don")) # 불용어 제거 ("english"는  a, an, the 등 관사, 전치사, 조사, 접속사, 불용어)
  # cps.tm <- tm_map(cps.tm, stemDocument) # Text stemming(어근만 추출한다) SnowballC패키지 설치
  return(cps.tm)
}

#---------

comments = readLines("./산출물_중간단계/COVID1.csv") #readLines : 한줄씩 읽어 문자열 벡터로 저장.
comments = readLines("./산출물_중간단계/COVID2.csv") #readLines : 한줄씩 읽어 문자열 벡터로 저장.
comments = readLines("COVID_19_children.csv") #readLines : 한줄씩 읽어 문자열 벡터로 저장.


comments <- str_replace_all(comments, "\\<[^>]+\\>"," ") #이모티콘 제거
comments <- str_replace_all(comments, "\\W"," ") #특수문자 제거 
comments <- str_replace_all(comments, "[가-힣]"," ") #한글 제거 
comments <- str_replace_all(comments, "schools","school") 
comments <- str_replace_all(comments, "masks","mask") 
comments <- str_replace_all(comments, "cases","case") 
comments <- str_replace_all(comments, "years","year") 
head(comments)

#-----------  앞에꺼 
tdm <- TermDocumentMatrix(clean.cps(comments))
inspect(tdm)
# dtm <- DocumentTermMatrix(clean.cps(comments))
# inspect(dtm)
m <- as.matrix(tdm) 
v <- sort(rowSums(m), decreasing = TRUE)

d_Rcran1 <- data.frame(word = names(v), freq = v)
rownames(d_Rcran1) <- 1:nrow(d_Rcran1)
head(d_Rcran1, 30)
#--------------  뒤에꺼 
d_Rcran2<- data.frame(word = names(v), freq = v)
rownames(d_Rcran2) <- 1:nrow(d_Rcran2)
head(d_Rcran2, 30)

d_all <- rbind(d_Rcran1, d_Rcran2)

d_all <- d_all %>% group_by(word) %>% 
  summarise(freq = sum(freq)) %>% 
  as.data.frame() %>% 
  arrange(desc(freq))
head(d_all,30)


word.plot <- wordcloud2(data = d_all)

#html형식으로 워드클라우드 파일 저장하기
saveWidget(word.plot, "./result/tmp3.html", selfcontained = F)

#tmp.html 파일을 png로 변환하기
webshot("./result/tmp3.html", "./result/wordcloud3.png", delay =20, vwidth = 1200, vheight=900)


# save(d_Rcran1, d_Rcran2, d_all, word.plot, file = "./R결과물/u.covid2.r")
#load("./R결과물/u.covid2.r")

save(d_Rcran1, d_Rcran2, d_all, word.plot, tdm, file = "./R결과물/u.covid2.r")
load("./R결과물/u.covid2.r")

save(d_all, file = "./R결과물/stem.r")
load("./R결과물/stem.r")


########################### 4. children + disaster safety

## 댓글 쪼개기 
disaster <- read.csv("children_disaster_safety.csv")
disaster1 <- disaster[1:25000,]
write.csv(disaster1, "disaster1.csv")
disaster2 <- disaster[25001:50000,]
write.csv(disaster2, "disaster2.csv")
disaster3 <- disaster[50001:79419,]
write.csv(disaster3, "disaster3.csv")

clean.cps <- function(comment){
  comments <- enc2utf8(comments)  # UTF-8로 변환
  cps <- Corpus(VectorSource(comments))
  cps.tm <- tm_map(cps, removePunctuation) # 구두점 제거
  cps.tm <- tm_map(cps.tm, removeNumbers) # 숫자제거
  cps.tm <- tm_map(cps.tm, tolower) #알파벳 모두 소문자로 변경 
  cps.tm <- tm_map(cps.tm, stripWhitespace) # 공백 제거 
  cps.tm <- tm_map(cps.tm,removeWords, c(stopwords("english"), stopwords("SMART"), 
                                         "don", "didn","lol")) # 불용어 제거 ("english"는  a, an, the 등 관사, 전치사, 조사, 접속사, 불용어)
  # cps.tm <- tm_map(cps.tm, stemDocument) # Text stemming(어근만 추출한다) SnowballC패키지 설치
  return(cps.tm)
}
#------------
comments = readLines("disaster1.csv") #readLines : 한줄씩 읽어 문자열 벡터로 저장.
comments = readLines("disaster2.csv")
comments = readLines("disaster3.csv")
comments = readLines("children_disaster_safety.csv")


comments <- str_replace_all(comments, "\\<[^>]+\\>"," ") #이모티콘 제거
comments <- str_replace_all(comments, "\\W"," ") #특수문자 제거 
comments <- str_replace_all(comments, "[가-힣]"," ") #한글 제거 
comments <- str_replace_all(comments, "schools","school") 
comments <- str_replace_all(comments, "masks","mask") 
comments <- str_replace_all(comments, "cases","case") 
comments <- str_replace_all(comments, "years","year") 

head(comments)

tdm <- TermDocumentMatrix(clean.cps(comments))
inspect(tdm)
m <- as.matrix(tdm) 
v <- sort(rowSums(m), decreasing = TRUE)

d_Rcran1 <- data.frame(word = names(v), freq = v)
rownames(d_Rcran1) <- 1:nrow(d_Rcran1)
head(d_Rcran1, 30)
#--------------
d_Rcran2<- data.frame(word = names(v), freq = v)
rownames(d_Rcran2) <- 1:nrow(d_Rcran2)
head(d_Rcran2, 30)
#--------------
d_Rcran3<- data.frame(word = names(v), freq = v)
rownames(d_Rcran3) <- 1:nrow(d_Rcran3)
head(d_Rcran3, 30)

d_all <- rbind(d_Rcran1, d_Rcran2, d_Rcran3)

d_all <- d_all %>% group_by(word) %>% 
  summarise(freq = sum(freq)) %>% 
  as.data.frame() %>% 
  arrange(desc(freq))
head(d_all,30)
dim(d_all)

word.plot <- wordcloud2(data = d_all)
word.plot

#html형식으로 워드클라우드 파일 저장하기
saveWidget(word.plot, "./result/tmp4.html", selfcontained = F)

#tmp.html 파일을 png로 변환하기
webshot("./result/tmp4.html", "./result/wordcloud4.png", delay =20,vwidth = 1200, vheight=900)

save(d_Rcran1, d_Rcran2, d_Rcran3, d_all, word.plot, file = "./R결과물/u.disaster2.r")
load("./R결과물/u.disaster2.r")

save(d_Rcran1, d_Rcran2, d_Rcran3, d_all, word.plot, tdm,  file = "./R결과물/u.disaster2.r")
load("./R결과물/u.disaster2.r")

#########
# 0128보고서 워드 클라우드 수정사항
setwd("C:/Users/hyerin/Google 드라이브(myfriands10@gmail.com)/Mixture_lab/youtube_산출물")
load("./R결과물/k.disaster.r")
library(wordcloud2)
a <- wordcloud2(word$two_top_100)
b <- wordcloud2(word$one_top_100)
c <- wordcloud2(word$top_100)

load("./R결과물/k.covid.r")
d <- wordcloud2(word$two_top_100)
e <- wordcloud2(word$one_top_100)
f <- wordcloud2(word$top_100)


#html형식으로 워드클라우드 파일 저장하기
saveWidget(f, "./k_cov_0128.html", selfcontained = F)

#tmp.html 파일을 png로 변환하기
webshot("./k_cov_0128.html", "./kcov2_0128.png", delay =20,vwidth = 1200, vheight=900)
webshot("./k_cov_0128.html", "./kcov1_0128.png", delay =20,vwidth = 1200, vheight=900)
webshot("./k_cov_0128.html", "./kcovall_0128.png", delay =20,vwidth = 1200, vheight=900)



################################################################
######   LDA   ###################################################
########################## topic_num
# 참고 사이트: https://rpubs.com/MNidhi/NumberoftopicsLDA

# Metrics used for Comparison
# Arun2010: The measure is computed in terms of symmetric KL-Divergence of salient distributions that are derived from these matrix factor and is observed that the divergence values are higher for non-optimal number of topics (maximize)
# 
# CaoJuan2009: method of adaptively selecting the best LDA model based on density.(minimize)
# 
# Griffths: To evaluate the consequences of changing the number of topics T, used the Gibbs sampling algorithm to obtain samples from the posterior distribution over z at several choices of T(minimize)
#    -> 논문. Thomas L. Griffiths and Mark Steyvers. 2004. Finding scientific topics. Proceedings of the National Academy of Sciences 101, suppl 1: 5228–5235. http://doi.org/10.1073/pnas.0307752101
# Perplexity:The most common way to evaluate a probabilistic model is to measure the log-likelihood of a held-out test set; Perplexity is a measurement of how well a probability distribution or probability model predicts a sample


#### number of topics 
## 1. for korean keyword

topic_num <- function(nouns){
  
  for(i in 1:length(nouns)){
    nouns[[i]] <- c(nouns[[i]][nchar(nouns[[i]])>=2])
  }
  
  nouns_1 <- list()
  for(i in 1:length(nouns)){
    nouns_1[[i]] <- paste(nouns[[i]], collapse = " ")
  }
  
  z <- Corpus(VectorSource(nouns_1)) 
  zz <- DocumentTermMatrix(z)
  rowTotals <- apply(zz , 1, sum)
  zz <- zz[rowTotals> 0, ]
  
  ##plot the metrics to get number of topics
  result <- FindTopicsNumber(
    zz,
    topics = seq(from = 2, to = 20, by = 1),
    metrics = c("Griffiths2004", "CaoJuan2009"),
    method = "Gibbs",
    control = list(seed = 1234),
    mc.cores = 4L,
    verbose = TRUE
  )
  
  FindTopicsNumber_plot(result)
}

load("./R결과물/k.covid2.r")
load("./R결과물/k.disaster2.r")

topic_num(nouns)

## 2. for english keyword

topic_num_en <- function(tdm){
  
  word.count = as.array(rollup(tdm,2))   #매트릭스 행별 합계구하기
  word.order = order(word.count, decreasing = T)[1:1000] #많이 쓰인 단어 순서정리(단어번호)
  freq.word = word.order[1:1000]  #상위 1000개 단어만 재할당(단어번호)
  row.names(tdm[freq.word,])      #해당단어번호 단어 확인
  
  #### LDA에 input data 만들기 (DTM만든후 다시 dtm2ldaformat함수로 변환)
  dtm = as.DocumentTermMatrix(tdm[freq.word,])  #상위 1000개 단어만 뽑아서 DTM으로 변환
  zz <- dtm
  rowTotals <- apply(zz , 1, sum)
  zz <- zz[rowTotals> 0, ]
  
  ##plot the metrics to get number of topics
  result <- FindTopicsNumber(
    zz,
    topics = seq(from = 2, to = 40, by = 1),
    metrics = c("Griffiths2004", "CaoJuan2009"),
    method = "Gibbs",
    control = list(seed = 1234),
    mc.cores = 4L,
    verbose = TRUE
  )
  
  FindTopicsNumber_plot(result)
}

load("./R결과물/u.covid2.r")
load("./R결과물/u.disaster2.r")

topic_num_en(tdm)


########################## LDA
##stopword추가 : 때문, 하세 

myLDA <- function(nouns, stopword ,K){
  
  # for(i in 1:length(nouns)){
  #   nouns[[i]] <- str_trim(nouns[[i]])
  #   nouns[[i]] <- nouns[[i]][which(!nouns[[i]] %in% stopwords)]
  #   nouns[[i]] <- c(nouns[[i]][nchar(nouns[[i]])>=2])
  # }
  
  temp <- sapply(nouns, str_trim)
  temp1 <- sapply(temp, FUN = function(x) x[which(x != "")])
  temp2 <- sapply(temp1, FUN = function(x) x[which(!x %in% stopwords)])
  nouns <- sapply(temp2, FUN = function(x) c(x[nchar(x) >= 2]))
  
  
  ##### LDA #####
  
  # nouns_1 <- list()
  # for(i in 1:length(nouns)){
  #   nouns_1[[i]] <- paste(nouns[[i]], collapse = " ")
  # }
  nouns_1 <- lapply(nouns, FUN = function(x) paste(x, collapse = " "))
  
  corpusLDA <- lexicalize(nouns_1)  # lexicalize : Generate LDA Documents from Raw Text Description
  
  
  word.proportion.top100 <- list()
  iter=0
  
  K <- K
  G <- 1000
  alpha <- 1
  eta <- 0.1
  
  set.seed(357)
  ldaModel=lda.collapsed.gibbs.sampler(corpusLDA$documents,K=K,vocab=corpusLDA$vocab,burnin=9999,num.iterations=G,alpha=alpha,eta=eta)
  
  top.words <- top.topic.words(ldaModel$topics, 20) # 각 주제의 상위 20개 단어만 출력
  
  topic.proportion <- ldaModel$topic_sums/sum(ldaModel$topic_sums)
  
  # word.proportion <- matrix(ncol=K,nrow=ncol(ldaModel$topics))
  # for(j in 1:K){
  #   word.proportion[,j] <- ldaModel$topics[j,]/ldaModel$topic_sums[j] 
  # }
  word.proportion <- (ldaModel$topics * as.vector(1/ldaModel$topic_sums)) %>%  t
  
  row.names(word.proportion) <- attributes(ldaModel$topics)$dimnames[[2]]
  
  topic.pct <- paste0(round(topic.proportion[,1],3)*100,"%")
  
  for(j in 1:K){  # x[order(-x)] : x를 내림차순 정렬 / x[order(x)] : 오름차순 정렬
    word.proportion.top100[[j]] <- round(word.proportion[c(order(-word.proportion[,j])),][1:100,j]*100,2)
  }

  
  # LDA분석 결과 상위 15개 단어 
  # new_LDA_result <- list()
  # for(i in 1:K){
  #   a<-head(word.proportion.top100[[i]],15)
  #   new_LDA_result[[i]] <- paste0(names(a)," (",a,"%)")
  # }
  a <- lapply(word.proportion.top100, function(x) head(x,15))
  new_LDA_result <- lapply(a, function(x) paste0(names(x)," (",x,"%)"))
  
  
  return(list(top.words=top.words, topic.proportion=topic.proportion, topic.pct=topic.pct,
              word.proportion.top100=word.proportion.top100,
              new_LDA_result=new_LDA_result, ldamodel = ldaModel))
}


# 유아 + 재난안전
load("./R결과물/k.disaster2.r")

LDA_result <- myLDA(nouns,stopwords, 15)
LDA_result <- myLDA(nouns,stopwords, 5)
LDA_result <- myLDA(nouns,stopwords, 5)
LDA_result <- myLDA_result4


LDA_result$top.words
LDA_result$topic.proportion
LDA_result$topic.pct
LDA_result$word.proportion.top100
LDA_result$new_LDA_result

k.dis.list = list()
iter = 1
for(i in 3:15){
  k.dis.list[[iter]] <- myLDA(nouns, stopwords,i)
  iter = iter + 1
  print(iter)
}

save(k.dis.list, nouns, stopwords, file = "./R결과물/k.disaster.lda_list.r")

for(i in 3:15){
  print(i)
  print(LDA.result[[i]]$topic.pct)
}

load("./R결과물/k.disaster_lda.r")
save(comments, nouns, stopwords, word, LDA_result4, LDA_result5, LDA_result6,LDA.result, file = "./R결과물/k.disaster_lda.r")

# 코로나19 + 유아
load("./R결과물/k.covid2.r")
stopwords = c(stopwords,"송혜교","유아인","하세")
# LDA_result4 <- myLDA(nouns,stopwords, 4)
# LDA_result5 <- myLDA(nouns,stopwords, 5)
# LDA_result6 <- myLDA(nouns,stopwords, 6)
# LDA_result <- myLDA_result4


k.covid.list = list()
iter = 1
for(i in 3:20){
  k.covid.list[[iter]] <- myLDA(nouns, stopwords,i)
  iter = iter + 1
  print(iter)
}
save(k.covid.list, nouns, stopwords, file = "./R결과물/k.covid.lda_list.r")


load("./R결과물/k.covid_lda.r")
save(comments, nouns, stopwords, word, LDA_result4, LDA_result5, LDA_result6, LDA.result, file = "./R결과물/k.covid_lda.r")


###############################################
####### 02/04 키워드 3,4 LDA 다른 방법

myLDA2 <- function(tdm, stopword ,K){
  
  word.count = as.array(rollup(tdm,2))   #매트릭스 행별 합계구하기
  word.order = order(word.count, decreasing = T)[1:1000] #많이 쓰인 단어 순서정리(단어번호)
  freq.word = word.order[1:1000]  #상위 1000개 단어만 재할당(단어번호)
  # row.names(tdm[freq.word,])      #해당단어번호 단어 확인
  
  #### LDA에 input data 만들기 (DTM만든후 다시 dtm2ldaformat함수로 변환)
  dtm = as.DocumentTermMatrix(tdm[freq.word,])  #상위 1000개 단어만 뽑아서 DTM으로 변환
  
  
  ldaform = dtm2ldaformat(dtm, omit_empty=T)    #dtm을 LDA 돌릴수 있는 형식의 데이터로 변환
  
  ldaModel = lda.collapsed.gibbs.sampler(documents = ldaform$documents,
                                         K = K,
                                         vocab = ldaform$vocab,
                                         num.iterations = 5000, #posteria 업데이트 회수    
                                         burnin = 1000, #초기값 버닝 
                                         alpha = 0.01,  #파라미터 세팅 : 문서 내에서 토픽들의 확률분포 / 1 = 유니폼
                                         eta = 0.01 )   #파라미터 세팅 : 한 토픽내에 단어들의 확률분포 /
  

  
  # ldaModel=lda.collapsed.gibbs.sampler(corpusLDA$documents,K=K,vocab=corpusLDA$vocab,burnin=9999,num.iterations=G,alpha=alpha,eta=eta)
  
  top.words <- top.topic.words(ldaModel$topics, 20)
  
  topic.proportion <- ldaModel$topic_sums/sum(ldaModel$topic_sums)
  
  # word.proportion <- matrix(ncol=K,nrow=ncol(ldaModel$topics))
  # for(j in 1:K){
  #   word.proportion[,j] <- ldaModel$topics[j,]/ldaModel$topic_sums[j] 
  # }
  word.proportion <- (ldaModel$topics * as.vector(1/ldaModel$topic_sums)) %>%  t
  
  row.names(word.proportion) <- attributes(ldaModel$topics)$dimnames[[2]]
  
  topic.pct <- paste0(round(topic.proportion[,1],3)*100,"%")
  
  word.proportion.top100 <- list()
  # iter=0
  for(j in 1:K){
    word.proportion.top100[[j]] <- round(word.proportion[c(order(-word.proportion[,j])),][1:100,j]*100,2)
  }
  
  # LDA분석 결과 상위 15개 단어 
  # new_LDA_result <- list()
  # for(i in 1:K){
  #   a<-head(word.proportion.top100[[i]],15)
  #   new_LDA_result[[i]] <- paste0(names(a)," (",a,"%)")
  # }
  a <- lapply(word.proportion.top100, function(x) head(x,15))
  new_LDA_result <- lapply(a, function(x) paste0(names(x)," (",x,"%)"))
  
  return(list(top.words=top.words, topic.proportion=topic.proportion, topic.pct=topic.pct,
              word.proportion.top100=word.proportion.top100,
              new_LDA_result=new_LDA_result, dtm))
}


##### 3. COVID19+children
# 251 ~ 264줄 실행 후 
load("./R결과물/u.covid2.r")

stopwords = c(stopwords("english"),stopwords("SMART"),"don")
LDA_result4 <- myLDA2(tdm,stopwords, 4)
LDA_result5 <- myLDA2(tdm,stopwords, 5)
LDA_result6 <- myLDA2(tdm,stopwords, 6)
LDA_result <- LDA_result4

LDA_result$top.words
LDA_result$topic.proportion
LDA_result$topic.pct
LDA_result$word.proportion.top100
LDA_result$new_LDA_result
save(LDA_result4, LDA_result5, LDA_result6, file = "./R결과물/u.covid_LDA2.r")


u.covid.list = list()
iter = 1
for(i in 3:20){
  u.covid.list[[iter]] <- myLDA2(tdm, stopwords,i)
  iter = iter + 1
  print(iter)
}
save(u.covid.list, nouns, stopwords, file = "./R결과물/u.covid.lda_list.r")



######## 4. disaster + children
# 334 ~ 347줄 실행 후 코드 돌릴 것.
load("./R결과물/u.disaster2.r")
LDA_result4 <- myLDA2(tdm,stopwords, 4)
LDA_result5 <- myLDA2(tdm,stopwords, 5)
LDA_result6 <- myLDA2(tdm,stopwords, 6)
LDA_result <- LDA_result4

LDA_result$top.words
LDA_result$topic.proportion
LDA_result$topic.pct
LDA_result$word.proportion.top100
LDA_result$new_LDA_result

u.dis.list = list()
iter = 1
for(i in 3:20){
  u.dis.list[[iter]] <- myLDA2(tdm, stopwords,i)
  iter = iter + 1
  print(iter)
}
save(u.dis.list, nouns, stopwords, file = "./R결과물/u.dis.lda_list.r")

load("./R결과물/u.disaster_LDA.r")
save(LDA_result4, LDA_result5, LDA_result6,LDA.result, file = "./R결과물/u.disaster_LDA.r")


###########################################
################ LDA results (topic num : 3~20)

##### 1. 유아 + 재난안전
load("./R결과물/k.disaster2.r")
k.dis.list = list()
iter = 1
for(i in 3:20){
  k.dis.list[[iter]] <- myLDA(nouns, stopwords,i)
  iter = iter + 1
  print(iter)
}
save(k.dis.list, nouns, stopwords, file = "./R결과물/k.disaster.lda_list.r")

##### 2. 코로나19 + 유아
load("./R결과물/k.covid2.r")
k.covid.list = list()
iter = 1
for(i in 3:20){
  k.covid.list[[iter]] <- myLDA(nouns, stopwords,i)
  iter = iter + 1
  print(iter)
}
save(k.covid.list, nouns, stopwords, file = "./R결과물/k.covid.lda_list.r")

##### 3. COVID19+children
load("./R결과물/u.covid2.r")
u.covid.list = list()
iter = 1
for(i in 3:20){
  u.covid.list[[iter]] <- myLDA2(tdm, stopwords,i)
  iter = iter + 1
  print(iter)
}
save(u.covid.list, nouns, stopwords, file = "./R결과물/u.covid.lda_list.r")

######## 4. disaster + children
load("./R결과물/u.disaster2.r")
u.dis.list = list()
iter = 1
for(i in 3:20){
  u.dis.list[[iter]] <- myLDA2(tdm, stopwords,i)
  iter = iter + 1
  print(iter)
}
save(u.dis.list, nouns, stopwords, file = "./R결과물/u.dis.lda_list.r")


save(k.dis.list,k.covid.list, u.dis.list, u.covid.list, file = "./R결과물/lda_list4.r")


# topic.proportion according to the topic number
load("./R결과물/lda_list4.r")

lapply(k.dis.list, FUN = function(x) x$topic.proportion %>% t)
lapply(k.covid.list, FUN = function(x) x$topic.proportion %>% t)
lapply(u.dis.list, FUN = function(x) x$topic.proportion %>% t)
lapply(u.covid.list, FUN = function(x) x$topic.proportion %>% t)

































##########################################################
################## 이전


# COVID19+children
load("./R결과물/u.covid2.r")

# 참고 블로그 : http://textmining.kr/?p=432
# 빈출 단어만 간추리기  
word.count = as.array(rollup(tdm,2))   #매트릭스 행별 합계구하기
word.order = order(word.count, decreasing = T)[1:1000] #많이 쓰인 단어 순서정리(단어번호)
freq.word = word.order[1:1000]  #상위 1000개 단어만 재할당(단어번호)
row.names(tdm[freq.word,])      #해당단어번호 단어 확인

#### LDA에 input data 만들기 (DTM만든후 다시 dtm2ldaformat함수로 변환)
dtm = as.DocumentTermMatrix(tdm[freq.word,])  #상위 1000개 단어만 뽑아서 DTM으로 변환
dtm

ldaform = dtm2ldaformat(dtm, omit_empty=T)    #dtm을 LDA 돌릴있 있는 형식의 데이터로 변환

result.lda = lda.collapsed.gibbs.sampler(documents = ldaform$documents,
                                         K = 6,
                                         vocab = ldaform$vocab,
                                         num.iterations = 5000, #posteria 업데이트 회수    
                                         burnin = 1000, #초기값 버닝 
                                         alpha = 0.01,  #파라미터 세팅 : 문서 내에서 토픽들의 확률분포 / 1 = 유니폼
                                         eta = 0.01 )   #파라미터 세팅 : 한 토픽내에 단어들의 확률분포 /
attributes(result.lda)
dim(result.lda$topics)  #5개 토픽에 1000개의 단어 출현   
result.lda$topics
top.topic.words(result.lda$topics) #5개 각 토픽 별 상위 20개 단어
class(result.lda$topics)
dim(result.lda$document_sums)
result.lda$document_sums[,1]


# 토픽 주제별 %
topic.proportion <- result.lda$topic_sums/sum(result.lda$topic_sums)
topic.pct <- paste0(round(topic.proportion[,1],3)*100,"%")


K = 5
K = 4
K = 6
lda_cnt <- matrix(0, nrow = 20, ncol = 2*K)

for (i in 1:K){
  top.topic.words(result.lda$topics)[,i]
  words20 <- which(names(result.lda$topics[i,]) %in% top.topic.words(result.lda$topics)[,i])
  word_cnt <- result.lda$topics[i,words20]
  lda_cnt[,(2*i-1)] <- top.topic.words(result.lda$topics)[,i]
  lda_cnt[,2*i] <- word_cnt[top.topic.words(result.lda$topics)[,i]]
}
lda_cnt
# k=4
lda_cnt4 <- as.data.frame(lda_cnt)
lda_cnt4
# k=5
lda_cnt5 <- as.data.frame(lda_cnt)
lda_cnt5
# k=6
lda_cnt6 <- as.data.frame(lda_cnt)
lda_cnt6

# 단어별 빈도를 %로 표현
LDA_cnt6 <- lda_cnt6
for (i in 1:K) {
  LDA_cnt6[,2*i] <- round((as.numeric(levels(lda_cnt6[,2*i])[lda_cnt6[,2*i]]) / result.lda$topic_sums[i]) * 100, 2)
}
LDA_cnt6


LDA_CNT6 <- matrix(0, nrow = 20, ncol = K)
for (i in 1:K) {
  LDA_CNT6[,i] <- paste0(LDA_cnt6[,2*i-1]," (",LDA_cnt6[,2*i],"%)")
}

LDA_CNT4
LDA_CNT5
LDA_CNT6


save(comments, dtm, result.lda, lda_cnt4, lda_cnt5, lda_cnt6, 
     LDA_CNT4, LDA_CNT5, LDA_CNT6, file = "./R결과물/u.covid_lda.r")
load("./R결과물/u.covid_lda.r")


# disaster + children
load("./R결과물/u.disaster2.r")

word.count = as.array(rollup(tdm,2))   #매트릭스 행별 합계구하기
word.order = order(word.count, decreasing = T)[1:1000] #많이 쓰인 단어 순서정리(단어번호)
freq.word = word.order[1:1000]  #상위 1000개 단어만 재할당(단어번호)
row.names(tdm[freq.word,])      #해당단어번호 단어 확인

#### LDA에 input data 만들기 (DTM만든후 다시 dtm2ldaformat함수로 변환)
dtm = as.DocumentTermMatrix(tdm[freq.word,])  #상위 1000개 단어만 뽑아서 DTM으로 변환
dtm

ldaform = dtm2ldaformat(dtm, omit_empty=T)    #dtm을 LDA 돌릴있 있는 형식의 데이터로 변환

result.lda = lda.collapsed.gibbs.sampler(documents = ldaform$documents,
                                         K = 4,
                                         vocab = ldaform$vocab,
                                         num.iterations = 5000, #posteria 업데이트 회수    
                                         burnin = 1000, #초기값 버닝 
                                         alpha = 0.01,  #파라미터 세팅 : 문서 내에서 토픽들의 확률분포 / 1 = 유니폼
                                         eta = 0.01 )   #파라미터 세팅 : 한 토픽내에 단어들의 확률분포 /
attributes(result.lda)
dim(result.lda$topics)  #5개 토픽에 1000개의 단어 출현   
result.lda$topics
result.lda$topic_sums             #단어합계 : 전체토픽 중 어떤 토픽해당 당하는 단어가 제일많은가?

top.topic.words(result.lda$topics) #5개 각 토픽 별 상위 20개 단어
class(result.lda$topics)
dim(result.lda$document_sums)
result.lda$document_sums[,1]

# 단어별 빈도를 숫자(n)로 표현
K = 5
K = 4
K = 6
lda_cnt <- matrix(0, nrow = 20, ncol = 2*K)

for (i in 1:K){
  top.topic.words(result.lda$topics)[,i]
  words20 <- which(names(result.lda$topics[i,]) %in% top.topic.words(result.lda$topics)[,i])
  word_cnt <- result.lda$topics[i,words20]
  lda_cnt[,(2*i-1)] <- top.topic.words(result.lda$topics)[,i]
  lda_cnt[,2*i] <- word_cnt[top.topic.words(result.lda$topics)[,i]]
}
lda_cnt

# 단어별 빈도를 %로 표현
LDA_cnt4 <- lda_cnt4
for (i in 1:K) {
  LDA_cnt4[,2*i] <- round((as.numeric(levels(lda_cnt4[,2*i])[lda_cnt4[,2*i]]) / result.lda$topic_sums[i]) * 100, 2)
}
LDA_cnt4


LDA_CNT4 <- matrix(0, nrow = 20, ncol = K)
for (i in 1:K) {
  LDA_CNT4[,i] <- paste0(LDA_cnt4[,2*i-1]," (",LDA_cnt4[,2*i],"%)")
}

LDA_CNT4
LDA_CNT5
LDA_CNT6

# k=4
lda_cnt4 <- as.data.frame(lda_cnt)
lda_cnt4
# k=5
lda_cnt5 <- as.data.frame(lda_cnt)
lda_cnt5
# k=6
lda_cnt6 <- as.data.frame(lda_cnt)
lda_cnt6

save(comments, dtm, result.lda, lda_cnt4, lda_cnt5, lda_cnt6,
     LDA_CNT4, LDA_CNT5, LDA_CNT6 ,file = "./R결과물/u.disaster_lda.r")
load("./R결과물/u.disaster_lda.r")


























##########################################################################
# ##############################

# corpus(말뭉치) 전처리
clean.text <- function(text){
  text <- tolower(text)  # 대문자를 소문자로 
  text <- str_replace_all(text, "\\<[^>]+\\>","") #이모티콘 제거
  text <- str_replace_all(text, "\\W"," ") #특수문자 제거
  text <- str_replace_all(text, "schools","school") 
  text <- str_replace_all(text, "masks","mask") 
  text <- str_replace_all(text, "cases","case")
  text <- removeNumbers(text)
  text <- bracketX(text)
  text <- replace_number(text)
  text <- replace_abbreviation(text)
  text <- replace_contraction(text)
  text <- replace_symbol(text)
  text <- removePunctuation(text)
  text <- stripWhitespace(text)   #쓸데 없는 공란 제거
  text <- str_replace_all(text, "americans", "america")
  
  indexes <- which(text == "")
  if(length(indexes) > 0){
    text <- text[-indexes]
  } 
  return(text)
}




# corpus(말뭉치) 전처리
clean.text <- function(text){
  text <- tolower(text)  # 대문자를 소문자로 
  text <- str_replace_all(text, "\\<[^>]+\\>","") #이모티콘 제거
  text <- str_replace_all(text, "\\W"," ") #특수문자 제거
  text <- str_replace_all(text, "schools","school") 
  text <- str_replace_all(text, "masks","mask") 
  text <- str_replace_all(text, "cases","case")
  text <- replace_abbreviation(text) #약어를 원형으로 변경
  text <- replace_contraction(text)
  text <- replace_symbol(text)  #기호를 문자로 ex) @ becomes "at".
  text <- str_replace_all(text, "americans", "america")
  
  indexes <- which(text == "")
  if(length(indexes) > 0){
    text <- text[-indexes]
  } 
  return(text)
}


clean.cps <- function(comment){
  comments <- enc2utf8(comments)  # UTF-8로 변환
  cps <- Corpus(VectorSource(comments))
  cps.tm <- tm_map(cps, removePunctuation) # 구두점 제거
  cps.tm <- tm_map(cps.tm, removeNumbers) # 숫자제거
  cps.tm <- tm_map(cps.tm, tolower) #알파벳 모두 소문자로 변경 
  cps.tm <- tm_map(cps.tm, stripWhitespace) # 공백 제거 
  cps.tm <- tm_map(cps.tm,removeWords, c(stopwords("english"), stopwords("SMART"), "don")) # 불용어 제거 ("english"는  a, an, the 등 관사, 전치사, 조사, 접속사, 불용어)
  return(cps.tm)
}


# # 2. 데이터 전처리 -------------------------------------------------------
make_corpus <- function(text) {
  text_clean <- clean.text(text)
  # text_source <- VectorSource(text_clean)
  # text_corpus <- VCorpus(text_source)
  corpus <- clean.cps(text_clean)
  return(corpus)
}
# 
cps.all <- make_corpus(comments)

# 
# # 3. 말뭉치를 데이터프레임으로 변환 --------------------------------------
# 
word_freq <- function(corpus) {
  doc_tdm <- TermDocumentMatrix(corpus)
  doc_m <- as.matrix(doc_tdm)
  doc_term_freq <- sort(rowSums(doc_m), decreasing = TRUE)
  doc_word_freqs <- data.frame(Var1 = names(doc_term_freq),
                               Freq = doc_term_freq) 
  rownames(doc_word_freqs) <- 1:nrow(doc_word_freqs)
  return(doc_word_freqs)
}

word_freqs <- word_freq(cps.all)
head(word_freqs,30)

# word_freq <- function(corpus) {
#   doc_tdm <- TermDocumentMatrix(corpus)
#   doc_m <- as.matrix(doc_tdm)
#   doc_term_freq <- rowSums(doc_m)
#   doc_word_freqs <- data.frame(term = names(doc_term_freq),
#                                num = doc_term_freq) %>% arrange(desc(num))
#   return(doc_word_freqs)
# }


# 
# com.tx <- VectorSource(comments)
# com.tx <- VectorSource(com.tx)
# com.tx <- VCorpus(com.tx)
# com.tx <- clean.cps(com.tx)
# word_freqs <- word_freq(com.tx)
# 
# 
# 
# ########################## topic_num
# 
# # nouns : ???? ?????? ????
# 
# topic_num <- function(nouns){
#   
#   for(i in 1:length(nouns)){
#     nouns[[i]] <- c(nouns[[i]][nchar(nouns[[i]])>=2], nouns[[i]][nouns[[i]]=="??"|nouns[[i]]=="??"|nouns[[i]]=="??"|nouns[[i]]=="??"])
#   }
#   
#   nouns_1 <- list()
#   for(i in 1:length(nouns)){
#     nouns_1[[i]] <- paste(nouns[[i]], collapse = " ")
#   }
#   ##### ???? ???? ??占쏙옙 #####
#   
#   z <- Corpus(VectorSource(nouns_1)) 
#   zz <- DocumentTermMatrix(z)
#   rowTotals <- apply(zz , 1, sum)
#   zz <- zz[rowTotals> 0, ]
#   
#   result <- FindTopicsNumber(
#     zz,
#     topics = seq(from = 2, to = 15, by = 1),
#     metrics = c("Griffiths2004", "CaoJuan2009"),
#     method = "Gibbs",
#     control = list(seed = 1234),
#     mc.cores = 4L,
#     verbose = TRUE
#   )
#   
#   FindTopicsNumber_plot(result)
# }
# 
# topic_num(nouns)
# 
# 
# ########################## LDA
# 
# # nouns : ???? ?????? ????, K : ???? ??
# 
# LDA <- function(nouns, K){
#   
#   for(i in 1:length(nouns)){
#     nouns[[i]] <- c(nouns[[i]][nchar(nouns[[i]])>=2], nouns[[i]][nouns[[i]]=="??"|nouns[[i]]=="??"|nouns[[i]]=="??"|nouns[[i]]=="??"])
#   }
#   
#   ##### LDA #####
#   
#   nouns_1 <- list()
#   for(i in 1:length(nouns)){
#     nouns_1[[i]] <- paste(nouns[[i]], collapse = " ")
#   }
#   corpusLDA <- lexicalize(nouns_1)  
#   
#   ##### ???? ???? ?????玖? ?????? ???? #####
#   
#   word.proportion.top100 <- list()
#   iter=0
#   
#   K <- K
#   G <- 1000
#   alpha <- 1
#   eta <- 0.1
#   
#   set.seed(357)
#   ldaModel=lda.collapsed.gibbs.sampler(corpusLDA$documents,K=K,vocab=corpusLDA$vocab,burnin=9999,num.iterations=G,alpha=alpha,eta=eta)
#   
#   top.words <- top.topic.words(ldaModel$topics, 20)
#   
#   # ???? ??占쏙옙
#   topic.proportion <- ldaModel$topic_sums/sum(ldaModel$topic_sums)
#   
#   word.proportion <- matrix(ncol=K,nrow=ncol(ldaModel$topics))
#   for(j in 1:K){
#     word.proportion[,j] <- ldaModel$topics[j,]/ldaModel$topic_sums[j] 
#   }
#   row.names(word.proportion) <- attributes(ldaModel$topics)$dimnames[[2]]
#   
#   # ???? ?? ?輧? ??占쏙옙 top100
#   
#   for(j in 1:K){
#     word.proportion.top100[[j]] <- round(word.proportion[c(order(-word.proportion[,j])),][1:100,j]*100,2)
#   }
#   return(list(top.words=top.words, topic.proportion=topic.proportion, word.proportion.top100=word.proportion.top100))
# }
# 
# LDA_result <- LDA(nouns, 5)
# 
# ##### ?娩? ???활????? ???? ?????? #####
# 
# LDA_result$top.words
# 
# LDA_result$topic.proportion
# 
# LDA_result$word.proportion.top100
# word
# 
# 
# #-------------------------------------
# # 명사 추출
# ko.words <- function(doc){
#   d <- as.character(doc)
#   extractNoun(d)
# }
# 
# # 명사/ 형용사/ 동사 추출
# ko.words <- function(doc){
#   d <- as.character(doc)
#   pos <- paste(SimplePos09(d))
#   extracted <- str_match(pos, '([가-힣]+)/[NP]')
#   keyword <- extracted[,2]
#   keyword[!is.na(keyword)]
# }
# comment1 <- ko.words(comments)
# 
# # 말뭉치 만들기 : 단어가 몇번 쓰였는지 카운트 하는 TermDocumentMetrix 작성 
# cps <- VCorpus(VectorSource(comment1))
# 
# # TermDocumentMetrix 만들기
# tdm <- TermDocumentMatrix(cps,
#                           control = list(tokenize=ko.words,
#                                          removePunctuation = T,
#                                          removeNumbers = T,
#                                          stopwords = c('네','넹','넴','으',
#                                                        '음','움','오','헐','앜')))
# # tokenize : 형태소 분석
# #removePunctuation : 마침표 제거
# # removeNumbers : 숫자 제거
# # stopwords : 문자 제거
# tdm <- as.matrix(tdm)
# View(tdm)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
