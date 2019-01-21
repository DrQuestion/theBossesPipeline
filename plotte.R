args = commandArgs(trailingOnly=TRUE)
results=read.table(args[1], sep="\t", header=T, row.names = 1)
experiment_name=args[2]
if (length(args)==2) {
  args[3] = "vulcano.html"
}

results <- results[!is.na(results$padj),]

library(plotly)
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(rainbow(2))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  
magnitudes=log2(results$baseMean)
c=cRamp(magnitudes)
p=plot_ly(x=results$log2FoldChange,y = -log10(results$padj),
          text = ~paste('SYMBOL: ',row.names(results), '<br> Adjusted P Value: ', results$padj,'<br> log Fold Change: ', results$log2FoldChange, '<br> Base Mean: ', results$baseMean),hoverinfo = 'text',type = 'scatter',mode='markers',marker = list(color=c, opacity=0.5))%>%
  layout(title=experiment_name,xaxis=list(title='log2_Fold_Change'), yaxis=list(title='pAdjusted_P_Value'))

legend.plot <- plot_ly(hoverinfo='none') %>% 
  add_markers(x = 1, 
              y = seq_len(length(magnitudes)),
              showlegend = F, 
              marker = list(color=cRamp(1:length(magnitudes)))) %>%
  layout(
    annotations = list(
      list(x = 0.25, 
           y = 1, 
           text = "Expression Magnitude", 
           showarrow = F, 
           xref='paper', 
           yref='paper')),
    xaxis = list(
      title='',
      zeroline=F,
      showline=F,
      showticklabels=F,
      showgrid=F),
    yaxis=list(
      showgrid=F,
      showline=F,
      zeroline=F,
      tickmode='array',
      tickvals=c(length(magnitudes)*0.08,length(magnitudes)*0.5,length(magnitudes)*0.93),
      ticktext = c('LOW','MEDIUM','HIGH')))

p <- subplot(legend.plot, p, widths = c(0.25, 0.75), 
        titleX=TRUE, titleY=TRUE)
htmlwidgets::saveWidget(as_widget(p), args[3])