## Write tables to Excel
WriteToExcel <- function(tbls, fname) {
  require(openxlsx);
  
  wb <- createWorkbook();
  modifyBaseFont(wb, fontSize=14, fontName='Courier');
  
  for (i in 1:length(tbls)) {
    if (is.null(rownames(tbls[[i]]))) x <- data.frame(tbls[[i]], stringsAsFactors = FALSE, check.names = FALSE) else {
      x <- data.frame(ID=rownames(tbls[[i]]), tbls[[i]], stringsAsFactors = FALSE, check.names = FALSE); 
      if (any(duplicated(colnames(x)))) x <- data.frame(x, check.names = TRUE);
    }
      
    w <- sapply(1:ncol(x), function(j) max(nchar(colnames(x)[j]), median(nchar(as.character(x[[j]])))));
    addWorksheet(wb, sheetName = names(tbls)[i]);
    freezePane(wb, sheet = i, firstRow = TRUE, firstCol = TRUE);
    setColWidths(wb, sheet = i, cols = 1:ncol(x), widths = 2 + w);
    writeDataTable(wb, sheet = i, x = x, colNames = TRUE, rowNames = FALSE);
  }
  
  saveWorkbook(wb, file = fname, overwrite = TRUE);
}
