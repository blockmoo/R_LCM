library(rvest)
library(R6)
library(rJava)
library(hash)
library(Sejong)
library(KoNLP)
library(data.table)
library(reshape2)
library(reshape)
library(dplyr)
library(devtools)
library(wordcloud2)
library(RColorBrewer)
library(N2H4)
library(tm)
library(stringr)
library(tidyverse)
library(wordcloud)


#### KoNLP설치 : https://minimax95.tistory.com/entry/KoNLP-%ED%8C%A8%ED%82%A4%EC%A7%80-%EC%84%A4%EC%B9%98-%EC%98%A4%EB%A5%98-%ED%95%B4%EA%B2%B0-%EB%B0%A9%EB%B2%95
useSejongDic()

KoNLP::buildDictionary(
  ext_dic = c('sejong','woorimalsam')
)

### review 추출 함수
rev <- function(reviews) {
  reviews <- gsub("\\d+","",reviews) 
  reviews <- gsub("\\n","",reviews) 
  reviews <- substr(reviews,regexpr("점 중",reviews), regexpr("신고",reviews))
  reviews <- gsub("점 중\\t*","",reviews)
  reviews <- gsub("\\t*신","",reviews)
}

####### 댓글 추출 ---

url_base <- 'https://movie.naver.com/movie/point/af/list.nhn?st=mcode&sword=190010&target=&page='
all.reviews <- c()

for(page in 1:200){
  url <- paste(url_base,page,sep='')
  # read_html 함수를 사용하여 html 페이지를 htxt 변수에 저장
  htxt <- read_html(url,encoding="CP949")
  # html_nodes 함수를 사용하여 list_netizen class에 해당하는 부분을 table 변수에 저자
  table <- html_nodes(htxt,'.list_netizen') # html_nodes 함수를 사용하여 title class에 해당하는 부분을 content 변수에 저장 
  content <- html_nodes(table, 'td.title') %>% html_nodes("a.movie") %>% html_text()
   # html_text 함수를 사용하여 text 부분을 reviews 변수에 저장
   score <- html_nodes(table, "div.list_netizen_score") %>% html_nodes("em") %>% html_text()
   reviews <- html_nodes(table, "td.title")  %>% html_text()
   reviews <- rev(reviews)
   
  if(length(reviews)==0){break}
  all.reviews <- c(all.reviews, reviews)
  print(page)
}

head(all.reviews,10)
length(all.reviews)



#------------------------------
t1 <- write.csv(all.reviews,'naver.csv', row.names = F)
t2 <- table(t1)
t3 <- head(sort(t2,decreasing = T),30)
t3

tgt1 <- readLines("naver.csv", encoding = "euc-kr")
tgt1 <- tgt1[-1]
length(tgt1)  #2000
head(tgt1)

#### 형용사, 명사 추출 --------------------
# https://lightblog.tistory.com/55?category=733817

words <- function(doc){
  doc <- as.character(doc)
  doc2 <- paste(SimplePos22(doc))
  doc3 <- str_match(doc2, "([가-힣]+)/(PA)")  #NC|NQ|PA|PV
  doc4 <- doc3[,2]
  doc4[!is.na(doc4)]
}

word.list <- NULL
for (i in 1:length(tgt1)){
  word.list <- c(word.list, words(tgt1[i]))
}

word.uni <- unique(word.list)
View(word.uni)

two.word <- Filter(function(x){nchar(x)>=2},word.uni) 
two.word


###### 국립국어원DB 형용사
adj <- read.csv("./국립국어원/adj_word.csv")
head(adj)




######################################################################
##################### 
tgt1 <- readLines("naver.csv", encoding = "euc-kr")
tgt1 <- tgt1[-1]
comments <- str_replace_all(tgt1, "\\<.+\\>","") #이모티콘 제거
comments <- str_replace_all(comments, "\\W"," ") #특수문자 제거
comments <- str_replace_all(comments, "\\d"," ") #숫자 제거
head(comments,4)


#명사 추출 : 오류나는 문장들은 csv에서 수정 
nouns <- extractNoun(comments) 
length(nouns) #2000
head(nouns,4)


# 특수문자 삭제  
for(i in 1:length(nouns)){
  nouns[[i]] <- gsub('[~!@#$%&*()_+=?<>]','',nouns[[i]])
  nouns[[i]] <- gsub("\\[","",nouns[[i]])
  nouns[[i]] <- gsub('[ㄱ-ㅎ]','',nouns[[i]]) # ex) ㅋㅋㅋ,ㅎㅎㅎ
  nouns[[i]] <- gsub('(ㅜ|ㅠ)','',nouns[[i]])  # ex) ㅠㅠ, ㅜㅜ
  nouns[[i]] <- gsub("\\^","",nouns[[i]])
  nouns[[i]] <- gsub("[A-Z]","",nouns[[i]])
  nouns[[i]] <- gsub("[a-z]","",nouns[[i]])
}
head(nouns,4)
# View(nouns)

########################## word
### wordcloud
word_cloud <- function(nouns){
  
  wordcount <- table(unlist(nouns))
  wordcount <- sort(wordcount, decreasing = TRUE)
  
  df_word <- as.data.frame(wordcount, stringsAsFactors = F)
  df_word$Var1=as.character(df_word$Var1)

  two <- filter(df_word, nchar(Var1) >= 2) %>% # 두글자 이상만
            filter(Var1 != "  ")
  one <- filter(df_word, nchar(Var1) == 1)  # 한글자 단어만
  df_word2 <- rbind(two,one)
  
  all.word <- df_word2 %>%
    arrange(desc(Freq))
  one.word <- one %>%  
    arrange(desc(Freq))
  two.word <- two %>%
    arrange(desc(Freq))
  top.100 <- two.word %>% head(150)
  
  return(list(nouns=nouns, one.word = one.word,
              two.word=two.word, all.word=all.word, top.100 = top.100))
}

word <- word_cloud(nouns)
word$wordcloud
word$one.word %>% filter(Freq >=50)  #한글자 
word$two.word %>% filter(Freq >=50)  #두글자이상
head(word$top.100)
wordcloud2(word$two.word %>% head(150), size=1.2)

head(two)



##### 불용어 제거 (메모장은 ANSI로 저장)
stop.txt <- readLines("movie_stopwords.txt", encoding = "euc-kr" ) 
stop.txt <- str_replace_all(stop.txt, "\\<.+\\>","") #이모티콘 제거
stop.txt <- str_replace_all(stop.txt, "\\W","") #특수문자 제거
stop.txt

two.rm <- two %>% filter(!Var1 %in% stop.txt)
head(two.rm,30)



# 불용어 처리
# SMART 사전 내 있는 단어들은 전부 배제합니다.
all <-  word$two.word
stop <- all$Var1[!all$Var1 %in% stop.txt]
word$two.word <- word$two.word[which(word$two.word$Var1 %in% stop),]
dim(all)
dim(word$two.word)
length(stop)










###  수정 필요 -----------------------------------------------

t2 <- gsub("\\d+","",tgt1[2]) 
gsub("[[:punct:]]", "", t2)

words2 <- function(doc){
  d <- gsub("\\d+","",doc) 
  d2 <- gsub("[[:punct:]]", "", d)
  return(d2)
}

data1 <- words2(tgt1[2:4])
data2 <- extractNoun(data1)
tdm <- TermDocumentMatrix(Corpus(VectorSource(data2)),
                          control = list(removeNumbers = T,
                                         removePunctuation = T,
                                         stopwords = T,
                                         wordLength = c(2, Inf)))
findFreqTerms(tdm, lowfreq = 10)


#빈출단어만 간추리기
library(slam)
word.count = as.array(rollup(tdm,2))
word.order = order(word.count, decreasing = T)[1:1000]
freq.word = word.order[1:1000]
row.names(tdm[freq.word,1])


corp <- Corpus(VectorSource(data2))
corp <- tm_map(corp, removePunctuation)
corp <- tm_map(corp, removeNumbers)
corp <- tm_map(corp, tolower)

tdm <- TermDocumentMatrix(corp, control = list(wordLengths = c(2, Inf)))
# inspect frequent term
findFreqTerms(tdm)


as.matrix(tdm)   

options(mc.cores=1)

ko.words <- function(doc){
  d <- as.character(doc)
  extractNoun(d)
}
tdm <- TermDocumentMatrix(corp,
                          control=list(tokenize=ko.words,
                                       wordLengths=c(1,10)))
as.matrix(tdm)




##### 추후 확인인 ----------------------------
all.reviews <- tgt1
all.reviews <- all.reviews[!str_detect(all.reviews,"평점")]   # 수집에 불필요한 단어가 포함된 내용 제거
Encoding(all.reviews)
options(encoding="utf-8")

## 명사/형용사 추출 함수 생성

ko.words <- function(doc){
  d <- as.character(doc)
  pos <- paste(SimplePos09(d))
  extracted <- str_match(pos, '([가-힣]+)/[NP]')
  keyword <- extracted[,2]
  keyword[!is.na(keyword)]
}

options(mc.cores=1)    # 단일 Core 만 활용하도록 변경 (옵션)
cps <- Corpus(VectorSource(word.list))  
tdm <- TermDocumentMatrix(cps,   
                          control=list(tokenize=ko.words,   ## token 분류시 활용할 함수명 지정
                                       removePunctuation=T,
                                       removeNumbers=T,
                                       wordLengths=c(2, 6),  
                                       weighting=weightBin
                          ))  

#최종결과 확인

dim(tdm)
tdm.matrix <- as.matrix(tdm)
Encoding(rownames(tdm.matrix)) <- "UTF-8"

word.count <- rowSums(tdm.matrix)  ##각 단어별 합계를 구함
word.order <- order(word.count, decreasing=T)  #다음으로 단어들을 쓰인 횟수에 따라 내림차순으로 정렬
freq.words <- tdm.matrix[word.order[1:20], ] #Term Document Matrix에서 자주 쓰인 단어 상위 20개에 해당하는 것만 추출
co.matrix <- freq.words %*% t(freq.words)  #행렬의 곱셈을 이용해 Term Document Matrix를 Co-occurence Matrix로 변경

par(family="Apple SD Gothic Neo")   ## mac option
qgraph(co.matrix,
       labels=rownames(co.matrix),   ##label 추가
       diag=F,                       ## 자신의 관계는 제거함
       layout='spring',              ##노드들의 위치를 spring으로 연결된 것 처럼 관련이 강하면 같이 붙어 있고 없으면 멀리 떨어지도록 표시됨
       edge.color='blue',
       vsize=log(diag(co.matrix))*2) ##diag는 matrix에서 대각선만 뽑는 것임. 즉 그 단어가 얼마나 나왔는지를 알 수 있음. vsize는 그 크기를 결정하는데 여기 인자값으로 단어가 나온 숫자를 넘겨주는 것임. log를 취한것은 단어 크기의 차이가 너무 커서 log를 통해서 그 차이를 좀 줄여준것임. 



########## kini'n creations
tgt1 <- readLines("naver.csv", encoding = "euc-kr")
tgt1 <- tgt1[-1]
head(tgt1)
mp <- SimplePos09(tgt1)

m_df <- mp %>% 
  melt %>% 
  as.tibble %>% 
  select(3,1)
head(m_df)

m_df1 <- m_df %>% 
  mutate(noun=str_match(value, '([가-힣]+)/N')[,2]) %>% 
  na.omit %>% 
  filter(str_length(noun)==1) %>% 
  count(noun, sort=TRUE)
word_m_df1 <-m_df1 %>% 
  filter(n>=20)
View(word_m_df1)


m_df2 <- m_df %>% 
  mutate(noun=str_match(value, '([가-힣]+)/N')[,2]) %>% 
  na.omit %>% 
  filter(str_length(noun)>=2) %>% 
  count(noun, sort=TRUE) %>% 
  filter(n>=20)
View(m_df2)









