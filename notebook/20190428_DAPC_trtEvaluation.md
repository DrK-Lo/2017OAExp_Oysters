# Examination of Discriminant function for Treatment

Initially we created a discriminant function for treatment using all of the data. This gave as a function that could only moderately discriminate by treatment, indicating that there was was sample variation occuring within treatment rathering between treatment, preventing clear discrimination. One reason for this could be variation between timepoints within a treatment is larger than variation among treatment leveles. 

![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/img/DAPC_wholeSampleTrt_PC10.png)

Next step was to examine perform a DAPC using only samples from the first timepoint (9) to create a discriminant funcion and than predict where along this discriminant fucntion the two treatment samples from the second timepoint fall:

![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/img/DAPC_Day9Trt_Day80Predict.png)

We see that a DAPC effectively discriminates between treatment levels using only samples from tp 9. However, when we predict the placement of tp 80 samples, the two treatment levels converge and overlap in the middle of discriminant function space. Indicating the the variation that was used to discriminate between treatment levels at tp 9 is not effective at discriminating between treatment after extended exposure. **ONE** possible explaination for this is that expression patterns are converging over time, becoming more similar to each other, but not necessarily reverting back to some original state (p_9 400). This could be the result of some other variable that is impacting expression and is also correlated with the variation that was used to create the discriminant function for treatment for the tp 9 samples

As a check I performed a DAPC with only the day points to check if it was capable of discriminating between treatments. We see that it can (keeping all parameter the selection the same, 10 PCs). Similarly I mapped th early timepoint onto the discriminant function and found that the discriminant function did not do a good job of discriminating between between treatment at day 9. 

![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/img/DAPC_Day9Trt_Day80Predict.png)
