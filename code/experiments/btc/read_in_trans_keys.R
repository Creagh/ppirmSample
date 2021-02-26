#########################################################
# SCRIPT: reads in the transaction keys for each 
# transaction ID
#
# CREATES: data frame of transaction keys
#########################################################

################################
# READ IN TRANSACTION KEYS     #
################################

# Start the clock!
ptm <- proc.time()
trans_key <- fread('/Users/cdb2/ownCloud/data/fabian_btc_network/transaction_key_list.csv', sep=",", header=FALSE,
									 col.names=c('transaction_ID', 'transaction_key'))
# Stop the clock
proc.time() - ptm
