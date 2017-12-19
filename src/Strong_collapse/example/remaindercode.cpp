//  void setRowZero(double rowIndx) // sets the rowIndx to zero implicitly and setting the respective non-zero columns for domination check
	// {
 //  		vertDomnIndicator[rowIndx] = true;  // Setting the vertex to be dominated.
 //  		for (sparseRowMatrix::InnerIterator itCol(sparseRowMxSimplices,rowIndx); itCol; ++itCol)  // Iterate over the non-zero columns
 //     		if(not simpDomnIndicator[itCol.index()] && not colInsertIndicator[itCol.index()])	  // Checking if the column is already dominated(set zero) or inserted	
 //      		{
 //      			columnIterator.push(itCol.index());
 //      			colInsertIndicator[itCol.index()] = true ;
 // 			} 		
 //  	}
 //  void setColumnZero(double colIndx)  
 //  {
 //    simpDomnIndicator[colIndx] = true;		// Setting the simplex to dominated.
 //    for (sparseMatrix::InnerIterator itRow(sparseMxSimplices,colIndx); itRow; ++itRow)  // Iterate over the non-zero columns
 //      if(not vertDomnIndicator[itRow.index()] && not rowInsertIndicator[itRow.index()]) // Checking if the row is already dominated(set zero) or inserted	
 //      {  
 //        rowIterator.push(itRow.index());
 //        rowInsertIndicator[itRow.index()] = true;
 //      }
 //  }
  

 
 //  doubleVector readRow(double rowIndx) // Returns list of non-zero Columns of the particular rowIndx
 //  {
 //    doubleVector columns; 
 //      for (sparseRowMatrix::InnerIterator itCol(sparseRowMxSimplices,rowIndx); itCol; ++itCol)  // Iterate over the non-zero columns
 //        if(not simpDomnIndicator[itCol.index()])
 //            columns.push_back(itCol.index()); // inner index, here it is equal to it.row()
    
 //      return columns;
 //  }

	// doubleVector readColumn(double colIndx) // Returns list of non-zero Rows of the particular colIndx
	// {
	// 	doubleVector rows; 
 //  		for (sparseMatrix::InnerIterator itRow(sparseMxSimplices,colIndx); itRow; ++itRow)  // Iterate over the non-zero columns
 //     		if(not vertDomnIndicator[itRow.index()])
 //     			rows.push_back(itRow.index());                                                  // inner index, here it is equal to it.row()
 
 //  		return rows;
	// }