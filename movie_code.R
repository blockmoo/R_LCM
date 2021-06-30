#+--------------------------------------------------------------------------------------------+#
# 제 목 : 네이버 영화 랭킹 페이지에서 영화 코드, 영화 제목, 평점을 2000개 가져온다.
# 작성자 : 서혜린 
# 작성일 : 2020-02-08
# 출력값 : 영화코드값, 영화명, 영화평점이 들어있는 파일을 생성
# 파일위치 : D:/myfriands10_google/Mixture_lab/naver_movie/영화댓글/movie_list_p400.txt
#+--------------------------------------------------------------------------------------------+#

library(XML)
library(stringr)

#크롤링 시작#
url_base <- "http://movie.naver.com/movie/sdb/rank/rmovie.nhn?sel=pnt&date=20210208&page=" #네이버 영화 랭킹 목록

all.codes <- c() #영화 코드 목록
all.titles <- c() #영화 이름 목록
all.points <- c() #영화 평점 목록

### '네이버 영화 랭킹 페이지'에서 1위부터 100위까지의 영화 목록을 수집
# 1-40페이지 영화 목록 수집 (40페이지가 최대, 총 2000개 영화)

iter = 0
for(page in 1:40){
  url <- paste(url_base, page, sep='')
  txt <- readLines(url, encoding="euc-kr")
  
  movie_info <- txt[which(str_detect(txt, "class=\"tit5\""))+1] #tit5클래스 아래 1줄 아래에는 영화 고유코드와 제목이 있다.
  points <- txt[which(str_detect(txt, "class=\"tit5\""))+7] #tit5클래스 아래 7줄 아래에는 평점이 있다.
  
  #titles #print
  #points #print
  
  codes <- substr(movie_info, 40, 50) #일부 코드를 선택
  codes <- gsub("\\D","",codes) #코드 중 숫자만 남기고 지우기
  titles <- gsub("<.+?>|\t", "", movie_info) # 텍스트만 남기고 코드 지우기 (이렇게하면 소스코드인식을 안하는듯)
  points <- gsub("<.+?>|\t", "", points) # 텍스트만 남기고 코드 지우기
  
  all.codes <-  c(all.codes, codes) #영화 코드값 저장
  all.titles <- c(all.titles, titles) #영화 이름 저장
  all.points <- c(all.points, points) #영화 평점 저장
  iter = iter + 1
  cat("iteration =" , iter, "\n")
}


#txt 파일로 출력
x <- cbind(all.codes, all.titles, all.points)
colnames(x) <- c("code", "movie_title", "point")

#x <- data.frame(code=c(all.codes), movie_title=c(all.titles),point=c(all.points))
dim(x)
head(x)

write.table(x, "./영화댓글/movie_list_p400.txt", row.names = F)


x <- as.data.frame(x)
class(x)
a <- match.arg(unique(dat$EPSD_NM), x$movie_title, several.ok = T)
head(a)  #1986개 일치 
























