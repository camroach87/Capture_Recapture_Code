TC_Chao_window_size_10.png
produced by running calcCR for window=10
mtrxCapt <- mtrxCaptD
window <- 10
qplot(y=aChao_endfix[,1]) + ylab("N") + xlab("Occasion") + geom_line(y=sChao$y) + ggtitle("Chao population estimates over secondary occasions for window size 10.") + theme_bw()
ggsave(file=file.path(plotDir,"TC_Chao_window_size_10.png"),width=10,height=6)

============================================================================================================

