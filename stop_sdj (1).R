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

######### 0. 전처리
options(mc.cores=1)
setwd("D:/myfriands10_google/Mixture_lab/naver_movie")

tgt1 <- readLines("./영화댓글/naver_tenet.csv", encoding = "euc-kr")
tgt1 <- readLines("./영화댓글/naver_inception.csv", encoding = "euc-kr")
tgt1 <- readLines("./영화댓글/naver_avengers.csv", encoding = "euc-kr")
tgt1 <- readLines("./영화댓글/naver_frozen.csv", encoding = "euc-kr")


tgt1 <- tgt1[-1]
comments <- str_replace_all(tgt1, "\\<[^>]+\\>","") #이모티콘 제거
comments <- str_replace_all(comments, "\\W"," ") #특수문자 제거
comments <- str_replace_all(comments, "\\d"," ") #숫자 제거
comments <- gsub('[~!@#$%&*()_+=?<>]','',comments)
comments <- gsub("\\[","",comments)
comments <- gsub('[ㄱ-ㅎ]','',comments) # ex) ㅋㅋㅋ,ㅎㅎㅎ
comments <- gsub('(ㅜ|ㅠ)','',comments)  # ex) ㅠㅠ, ㅜㅜ
comments <- gsub("\\^","",comments)
comments <- gsub("[A-Z]","",comments)
comments <- gsub("[a-z]","",comments)
head(comments,4)


##### 불용어 제거 (메모장은 ANSI로 저장)
stop.txt <- readLines("movie_stopwords.txt", encoding = "euc-kr" ) 
stop.txt <- str_replace_all(stop.txt, "\\<[^>]+\\>","") #이모티콘 제거
stop.txt <- str_replace_all(stop.txt, "\\W","") #특수문자 제거
stop.txt

#보통명사, 형용사 추출 함수
words <- function(doc){
  doc <- as.character(doc)
  doc2 <- paste(SimplePos22(doc))
  doc3 <- str_match(doc2, "([가-힣]+)/(PA)")
  doc4 <- doc3[,2]
  doc4[!is.na(doc4)]
}



#2.말뭉치변환
doc <- Corpus(VectorSource(comments)) 

doc <- TermDocumentMatrix(doc, #3.TermDocumentMatrix변환
control=list(                   #아래로옵션을나열
tokenize=words,                 #미리만들어둔함수(보통명사추출)로문장을자름
removeNumbers=T,                #숫자제거
removePunctuation=T,            #문장부호제거
wordLengths=c(1, 10),            #1~5음절로이루어진단어만추출
stopwords = stop.txt            #불용어 제거
))

doc <- as.matrix(doc) 
#4.Matrix로변환
doc <- rowSums(doc)
#5.행의합(총빈도)
doc <- doc[order(doc, decreasing=T)]
#빈도역순정렬
doc <- as.data.frame(doc)
doc <- data.frame(Var1 = rownames(doc), Freq = doc$doc, stringsAsFactors = F)
one <- filter(doc, nchar(Var1) == 1)
one
two <- filter(doc, nchar(Var1) >= 2)
two
head(two)

save(tgt1, doc, one, two, file = "./댓글추출중간/tenet.r")
save(tgt1, doc, one, two, file = "./댓글추출중간/inception.r")
save(tgt1, doc, one, two, file = "./댓글추출중간/avengers.r")
save(tgt1, doc, one, two, file = "./댓글추출중간/frozen.r")




###### 국립국어원DB 형용사
adj <- read.csv("./국립국어원/adj_word.csv")
adj$x <- as.character(adj$x)
# 형용사는 앞에 3글자만 따오기
adj.two <- str_sub(adj$x, 1,3)
adj.two.u <- unique(adj.two)
# 빈출어는 앞에 2글자만 따오기
two2 <- str_sub(two$Var1, 1,2)
two2.u <- unique(two2)

# 빈출어 2글자가 형용사 3글자에 매치했을 때 매칭되는 형용사만 가져옴.
x <- match.arg(two2, adj.two, several.ok = T)
length(x)

x2 <- match.arg(two2.u ,adj.two.u , several.ok = T)
length(x2)
head(x2,50)

# 빈출어 2글자가 형용사에 매치했을 때 매칭되는 형용사만 가져옴.
x3 <- match.arg(two2.u ,adj$x, several.ok = T)
length(x3)
head(x3,50)


# problem : 빈출어 2글자에 "하다" 가 포함되어있어서 형용사에 "하다"가 들어있으면 모두 출력
# -> "하다" 제거
which(two2 == "하다") #126
which(two2 == "이다") #466

two2 <- two2[-c(126,466)]

x2 <- match.arg(two2 ,adj.two.u , several.ok = T)
x2 <- unique(x2)
length(x2) #146개
head(x2,50)

which(two2.u == "하다")   #126
which(two2.u == "이다")   #466
two2.u <- two2.u[-c(126,466)]

x3 <- match.arg(two2.u ,adj$x, several.ok = T)
length(x3) #141개
head(x3,50)


##### two에서 Freq >= 10인 단어들만 추출
two.10 <- two %>% filter(Freq >= 10)
two10 <- str_sub(two.10$Var1, 1,2)
two10.u <- unique(two10)

which(two10.u == "하다")   #181
two10.u <- two10.u[-181]
which(two10.u == "이다")   #없음

x10 <- match.arg(two10.u ,adj$x, several.ok = T)
length(x10)  # 24개
x10
#problem: "다소", "사실"이라는 단어가 "다소곳하게", "사실무근하다"와 매칭됨.


###### 국립국어원DB 형용사
adj <- read.csv("./국립국어원/adj_word.csv")

###### 국립국어원DB 명사 
noun <- read.csv("./국립국어원/noun_word.csv")

###### 국립국어원DB 명사 
verb <- read.csv("./국립국어원/verb_word.csv")


capture <- function(wordtype) {
  wordtype$x <- as.character(wordtype$x)
  wordtype <- wordtype %>% filter(nchar(x) >= 2)
  
  # 앞에 3글자만 따오기
  word.two <- str_sub(wordtype$x, 1,4)
  word.two.u <- unique(word.two)
  
  x2 <- match.arg(two2.u ,word.two.u , several.ok = T)
  x3 <- match.arg(two2.u ,wordtype$x, several.ok = T)
  
  ##### two에서 Freq >= 10인 단어들만 추출
  x10 <- match.arg(two10.u ,wordtype$x, several.ok = T)
  return(list(word = wordtype, word.two = word.two, word.two.u = word.two.u,
              match2 = x2, match3 = x3, match10 = x10))
  
}

two10.u %>% head(30)


noun.cap <- capture(noun)
noun.cap %>%  sapply(head)
noun.cap %>% sapply(length)
noun.cap$match10


adj.cap <- capture(adj)
adj.cap %>% sapply(head)
adj.cap %>% sapply(length)
adj.cap$match10

verd.cap <- capture(verb)
verd.cap %>% sapply(head)
verd.cap %>% sapply(length)
verd.cap$match10



##############################################
###### 분위기 형용사 
##############################################
mood <- read.csv("./국립국어원/mood_word.csv")
head(mood)
names(mood) <- c("en.word", "word")


###### 단어들
load("./댓글추출중간/tenet.r")
load("./댓글추출중간/frozen.r")
load("./댓글추출중간/avengers.r")
load("./댓글추출중간/inception.r")



#freq >= 10인 단어들만 따오기
two.10 <- two %>% filter(Freq >= 10)
two10 <- str_sub(two.10$Var1, 1,2)
two10.u <- unique(two10)



# 빈출어는 앞에 2글자만 따오기
two2 <- str_sub(two$Var1, 1,2)
two2.u <- unique(two2)


capture2 <- function(wordtype) {
  wordtype$word <- as.character(wordtype$word)
  wordtype <- wordtype %>% filter(nchar(word) >= 2)
  
  # 앞에 3글자만 따오기
  word.two <- str_sub(wordtype$word, 1,3)
  word.two.u <- unique(word.two)
  
  x2 <- match.arg(two2.u ,word.two.u , several.ok = T)
  x2 <- match.arg(x2 ,wordtype$word , several.ok = T)
  
  x3 <- match.arg(two2.u ,wordtype$word, several.ok = T)
  x3 <- match.arg(x3 ,wordtype$word , several.ok = T)
  
  ##### two에서 Freq >= 10인 단어들만 추출
  x10 <- match.arg(two10.u ,word.two.u, several.ok = T)
  x10 <- match.arg(x10 ,wordtype$word , several.ok = T)
  return(list(word = wordtype, word.two = word.two, word.two.u = word.two.u,
              match2 = x2, match3 = x3, match10 = x10))
  
}


mood.cap <- capture2(mood)
mood.cap %>% sapply(head)
mood.cap %>% sapply(length)
mood.cap$match2

save(mood, two2, two, two2.u, capture2, mood.cap, file="./댓글추출중간/mood_tenet.r")
save(mood, two2, two, two2.u, capture2, mood.cap, file="./댓글추출중간/mood_frozen.r")
save(mood, two2, two, two2.u, capture2, mood.cap, file="./댓글추출중간/mood_avengers.r")
save(mood, two2, two, two2.u, capture2, mood.cap, file="./댓글추출중간/mood_inception.r")

# https://kerpect.tistory.com/132?category=876378


###################################
### 파이썬 파일 이용
###################################
setwd("D:/myfriands10_google/Mixture_lab/naver_movie/영화댓글")

fil <- read.csv("py_tenet.csv", col.names = F, encoding = "UTF-8")
names(fil) <- "Freq"

fil2 <- data.frame(word = rownames(fil), type = rownames(fil), Freq = fil, stringsAsFactors = F)
rownames(fil2) <- NULL
head(fil2)

for (i in 1:nrow(fil2)){
  fil2[i,1]  <- str_sub(fil2[i,1],3,(str_locate(fil2[i,1],",")[1]-2))
  fil2[i,2]  <- str_sub(fil2[i,2],(str_locate(fil2[i,2],",")[1]+3),-3)
}
head(fil2)






















