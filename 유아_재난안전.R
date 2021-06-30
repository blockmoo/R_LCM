# github : https://github.com/Leeyua-airim/R_AIRIM/blob/master/R_Selenium_Youtube/%EC%85%80%EB%A0%88%EB%8B%88%EC%9B%80_%EB%8B%A4%EC%9A%B4%EB%A1%9C%EB%93%9C_%EC%A3%BC%EC%86%8C
# 참고 유튜브 : https://youtu.be/aoCBKVPKJgk
#CMD 창 열어서 selenium 파일로 경로 변경 (cd 경로)
#java -Dwebdriver.gecko.driver="geckodriver.exe" -jar selenium-server-standalone-3.141.59.jar -port 4445

library(httr)
library(rvest)
# install.packages("RSelenium")
library(RSelenium)
library(dplyr)
library(stringr)



remD <- remoteDriver(port = 4445L, #포트번호 입력
                     browserName = "chrome") #사용할 브라우저

remD$open() #서버에 연결

title_youtube <- "코로나19+%2B+유아&sp=EgIIBQ%253D%253D" #검색어 & 필터 : 올해 , 동영상 
remD$navigate(paste0("https://www.youtube.com/results?search_query=",title_youtube)) #해당 홈페이지로 이동
# 고정된 주소 + 검색어를 paste()로 묶어줌

# 홈페이지 스크롤
remD$executeScript("window.scrollTo(0,500)")
for(i in 1:100){
  remD$executeScript("window.scrollTo(1000000,1500000)")
  Sys.sleep(1) #1초간 쉬게하는 코드 
}

html <- remD$getPageSource()[[1]] #페이지의 소스 가져오기   
html <- read_html(html) #페이지의 소스 읽어오기
html_url = html %>% html_nodes('#video-title')
href = html_attr(html_url,'href') #선택된 노드를 텍스트 화 
date <- html %>% html_nodes('#metadata-line') %>% 
  html_nodes('span') %>%   html_text()


url <- c()
for(i in 1:length(href)){
  url = c(url,paste0("https://www.youtube.com",href[i]))
}

length(url)  
View(url)

save(href, remD, title_youtube, html,html_url,url,file="usa_covid19_kids.r")
load("./youtube_산출물/korea_covid19_kids.r")

length(url)  #549개 
comments_result <- c()
iter = 0


for(j in 1:length(url)){
  remD$navigate(url[j])
  Sys.sleep(5)
  # 재생 및 일시정지를 위한 코드 
  btn <- remD$findElement(using="css selector", 
                          value = '.html5-main-video')
  btn$clickElement()
  Sys.sleep(5)
  
  remD$executeScript("window.scrollTo(2000,3000)") #(top:0, left:500) - top은 세로위치, left는 가로위치 
  Sys.sleep(8)
  html <- remD$getPageSource()[[1]] # 페이지 소스 가져오기 
  html <- read_html(html)
  
  # 댓글사용이 중지된 영상은 pass
  # ex) https://www.youtube.com/watch?v=dtjnrtmZiOU
  try(if(comment_num <- html %>% html_nodes("#message") %>% html_nodes("span") %>% html_text() == "댓글이 사용 중지되었습니다. "){
    iter = iter + 1
    cat("iteration =" , iter, "\n")
    next}, silent = T)
  
  # 댓글 갯수만큼만 for문 실행 
  comment_num <- html %>% html_nodes("#count") %>% html_nodes("yt-formatted-string") %>% html_text() 
  comment_num <- gsub(",","",str_sub(comment_num, 4,-2)) #댓글 1,125개 -> 1125
  scroll_num <- min(ceiling(as.integer(comment_num)/5),100)
  
  for(i in 1:scroll_num){
    remD$executeScript("window.scrollTo(1000000,1500000)")
    Sys.sleep(3)
    #new_height = remD$executeScript("return document.getElementById('hplogo');",args=list())
  }
  
  html <- remD$getPageSource()[[1]] # 페이지 소스 가져오기 
  html <- read_html(html)
  
  comments <- html %>% html_nodes("#content-text") %>% html_text() #선택된 노드를 텍스트화 
  comments <- gsub("\n", "", comments) 
  comments <- trimws(comments) #공백 제거 
  
  comments_result <- c(comments_result, comments)
  
  iter = iter + 1
  cat("iteration =" , iter, "\n")
}
# 
View(comments_result[1:10])
length(comments_result)
write.csv(comments_result, "./youtube_산출물/코로나19_유아.csv", row.names = F)










#-----------------------

iter95 <- read.csv("./youtube_산출물/COVID_19_children_iter11.csv")
iter96_166 <- read.csv("./youtube_산출물/COVID_19_children_iter252_570.csv")
iter167_559 <- read.csv("./youtube_산출물/COVID_19_children_iter12~250.csv")

dim(iter95) #4003개
dim(iter96_166) #1451개
dim(iter167_559) #11579개 
youtube_comments <- rbind(iter95, iter96_166, iter167_559) # 17033개 
write.csv(youtube_comments, "./youtube_산출물/COVID_19_children.csv", row.names = F)



### 영상 날짜 가져오는 코드
# url <- "https://www.youtube.com/watch?v=awqdBU7F3GU"
# remD$open()
# remD$navigate(url)
# html <- remD$getPageSource()[[1]]
# html <- read_html(html)
# text <- html %>% html_nodes('#date') %>% html_nodes('yt-formatted-string') %>% html_text()
# text

