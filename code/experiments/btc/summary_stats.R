################################################
# Compute SUMMARY STATS of the BTC user network#
################################################

library(ggplot2)

# STORAGE
object.size(df) # 1352812624 bytes ~= 1.35 GB

# NETWORK SIZE
nrow(df) # 37,450,461 rows
(num_tkeys <- length(unique(df$transaction_key))) # 15,898,622 unique transactions keys
# Note that multiple transactions (rows) can be associated with a single transaction key 

length(unique(df$from)) # 5,436,669 unique senders
length(unique(df$to)) # 6,335,769 unique recipients
(num_users <- length(unique(c(df$from, df$to)))) # 6,336,769 unique users IDs in the dataset

num_tkeys / num_users # 2.508948 transaction keys per unique user ID
nrow(df) / num_users # 5.910025 transactions (rows) per unique user ID

# AMOUNT SUMMARIES
summary(df$amount)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#    0.0      0.0      0.1     37.8      1.9  500000.0 
hist(df$amount)

# 1 + log scale
ggplot(df, aes(x=amount)) +
	stat_density(aes(y=..count..), color="black", fill="blue", alpha=0.3) +
	scale_x_continuous(trans="log1p", expand=c(0,0)) +
	theme_bw()
# extract the bottom 90% of transaction amounts
df2 <- df[df$amount < quantile(df$amount, 0.9), ]
hist(df2$amount)