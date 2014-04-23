TC_Chao_window_size_10.png
produced by running CR_RobustDesign for window=10
qplot(y=aChao_endfix[,1]) + ylab("N") + xlab("Occasion") + geom_line(y=sChao$y, colour="red") + ggtitle("Chao population estimates over secondary occasions for window size 10.") + theme_bw()

============================================================================================================

