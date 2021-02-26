library(twitteR)

consumer_key <- "4AY3bTm6lqHrFSZCbfj45bSAO"
consumer_secret <- "OIGRYj6bJ3OffmNzS6ZlORqxJNhT9cKVrJYugsWc9Z6lOCffzr"
access_token <- "2983304497-sKWAoVkU9lY42LIXDWEHaXC0SYQcejvvsOootUF"
access_secret <- "rWhfndcDH8fBPh3Oqa3TYTiaWmjWckSdLysPBAdnJuv7t"

setup_twitter_oauth(consumer_key, consumer_secret, access_token, access_secret)

ubc <- getUser("UBCStatistics")
friends.object <- lookupUsers(ubc$getFriendIDs())
followers.object <- lookupUsers(ubc$getFollowerIDs())

friends <- sapply(friends.object,name)
followers <- sapply(followers.object,name)

# Create a data frame that relates friends and followers to you for expression in the graph
relations <- merge(data.frame(User='UBCStatistics', Follower=friends),
									 data.frame(User=followers, Follower='UBCStatistics'), all=T)

n <- 50 # select a subset of followers to display
followers2 <- sapply(followers.object[1:n],name)
relations2 <- merge(data.frame(User='UBCStatistics', Follower=friends),
										data.frame(User=followers2, Follower='UBCStatistics'), all=T)




###########
tweets <- searchTwitter("ghomeshi", n=100, lang="en", since="2016-03-1")
# Transform tweets list into a data frame
tweets.df <- twListToDF(tweets)

###########

library("wordcloud")
library("tm")

#the cainfo parameter is necessary only on Windows
r_stats <- searchTwitter("micahgoldberg", n=1500)
#should get 1500
length(r_stats)
#[1] 1500

#save text
r_stats_text <- sapply(r_stats, function(x) x$getText())

#create corpus
r_stats_text_corpus <- Corpus(VectorSource(r_stats_text))

#clean up
r_stats_text_corpus <- tm_map(r_stats_text_corpus, content_transformer(tolower)) 
r_stats_text_corpus <- tm_map(r_stats_text_corpus, removePunctuation)
r_stats_text_corpus <- tm_map(r_stats_text_corpus, function(x)removeWords(x,stopwords()))
wordcloud(r_stats_text_corpus)

pal2 <- brewer.pal(8,"Dark2")
wordcloud(r_stats_text_corpus,min.freq=1,max.words=100, random.order=T, colors=pal2)
