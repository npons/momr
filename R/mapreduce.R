
#-------------------------------------------------------------------
# New parallel computing procedures. To be re-engeniered with an indexing
# approach and architecture of the files.
#-------------------------------------------------------------------

#' \code{splitData} 
#' @title splitData
#' @export
#' @description This function splits a dataframe object on a given number of equally sized shares
#' and saves them in the disk as RData objects with an incremental name.
#' @author Edi Prifti
#' @param data : the dataframe to be split
#' @param shares : the number of shares, default=30
#' @param pattern : the name of the split files in the disk preceding the incremental number, default="data_part"
#' @param rows : logical parameter, if default rows=="TRUE" the rows will be split
#' @return does not return anything
splitData <- function (data, shares = 30, pattern = "data_part", rows = TRUE) {
  if (rows) {
    genes_per_share <- ceiling(nrow(data)/shares)
    print(paste("There will be", genes_per_share,"elements per share."))
    coords <- as.data.frame(matrix(NA, ncol = 2, nrow = length(seq(1, nrow(data), genes_per_share))))
    coords[, 1] <- seq(1, nrow(data), genes_per_share)
    tmp <- seq(genes_per_share, nrow(data), genes_per_share)
    if(length(tmp) == nrow(coords)){
      coords[, 2] <- tmp
    } else {
      coords[, 2] <- c(tmp, nrow(data))
    }
    colnames(coords) <- c("start", "stop")
    for (i in 1:nrow(coords)) {
      print(paste("Cutting share", i))
      data_part <- data[coords$start[i]:coords$stop[i], ]
      saveRDS(data_part, file = paste(pattern, "_", i, ".rda", sep = ""), compress=FALSE)
    }
  }
  else {
    samples_per_share <- ceiling(ncol(data)/shares)
    coords <- as.data.frame(matrix(NA, ncol = 2, nrow = length(seq(1, ncol(data), samples_per_share))))
    coords[, 1] <- seq(1, ncol(data), samples_per_share)
    coords[, 2] <- c(seq(samples_per_share, ncol(data), samples_per_share), ncol(data))[1:nrow(coords)]  # modified elc 20150928
    if(coords[nrow(coords),1] == coords[nrow(coords),2]){ 
      coords <- coords[-nrow(coords), ] # if no elements at the last row delete
      coords[nrow(coords), 2] <- ncol(data) # close the edge
    }
    colnames(coords) <- c("start", "stop")
    for (i in 1:nrow(coords)) {
      print(paste("Cutting share", i))
      data_part <- data[, coords$start[i]:coords$stop[i]]
      saveRDS(data_part, file = paste(pattern, "_", i, ".rda", sep = ""), compress=FALSE)
    }
  }
}

#' \code{launchTask} 
#' @title launchTask
#' @export
#' @description This function distributes the calculations as separate processes in a multi-thread server.
#' @author Edi Prifti
#' @param input : a folder containing elements that are to be processed
#' @param output : a folder where the processed results will be writen
#' @param script : the R script that is to be executed in parallel
#' @param nb_proc : number of processes to run the task
#' @param resume : resume the computation for the remaining shares (default=TRUE)
#' @param pattern : the pattern of the files to be used in the // computation procedure
#' @return nothing
#' @note the number of processors may be specified. TODO: add the argument
launchTask <- function (input, output, script, nb_proc=NULL, resume=TRUE, pattern="data_part") {
  if(resume){# If we want to resume the process
    files.input <- dir(input, pattern = paste(pattern, sep = "."))
    files.output <- dir(output, pattern = paste(pattern, sep = "."))
    remaining <- files.input[!files.input %in% files.output]
    write.table(x=paste(input,remaining,sep="/"), file = "remaining.txt", row.names = F, quote = F, col.names = F)
    if(length(remaining)==length(files.input)){ 
      print(paste("Launching job",script)) 
    } else{
      print(paste("Resuming job",script,".", length(remaining),"shares remaining"))  
    }
    
    if(is.null(nb_proc)){
      # if nb_proc is NULL, all available processors will be used.
      command <- paste("cat remaining.txt | xargs -l -i --max-procs=`grep -c '^processor' /proc/cpuinfo` Rscript ", 
                       script, " {} ", output, "/ &>/dev/null &", sep = "")
    }else{
      command <- paste("cat remaining.txt | xargs -l -i --max-procs=",nb_proc," Rscript ", 
                       script, " {} ", output, "/ &>/dev/null &", sep = "")
    }
  }else{
    if(is.null(nb_proc)){
      # if nb_proc is NULL, all available processors will be used.
      command <- paste("nohup ls -1 ", input, "/* | xargs -l -i --max-procs=`grep -c '^processor' /proc/cpuinfo` Rscript ", 
                       script, " {} ", output, "/ &>/dev/null &", sep = "")
    }else{
      command <- paste("nohup ls -1 ", input, "/* | xargs -l -i --max-procs=",nb_proc," Rscript ", 
                       script, " {} ", output, "/ &>/dev/null &", sep = "")
    }
  }
  system(command) # launch command
  system("rm -fr remaining.txt") # clean space
}


#' \code{watchProgress} 
#' @title watchProgress
#' @export
#' @description This function prints the progress % of the overall threads that are treated. It is based
#' on the input and output folders by comparing the number of output objects.
#' @author Edi Prifti
#' @param input : a folder containing elements that are to be processed
#' @param output : a folder where the processed results will be writen
#' @return nothing
watchProgress <- function(input, output){
  progress <- 0
  while(progress!=100){
    input.files <- dir(input)
    output.files <- dir(output)
    flush.console()
    cat(progress <- length(output.files)/length(input.files)*100,"\r")
    flush.console()
    Sys.sleep(5)
    # TODO add chack stopping and wrap this into a function
  }
}


#' \code{mapResults} 
#' @title mapResults : merge data from // computing (the reduce process)
#' @export
#' @description This function loads the results once finished and merges them together in one dataframe. For the moment this is a slow procedure and we are testing ways to accelerate the binding.
#' @author Edi Prifti
#' @param folder : the folder where the results are found
#' @param pattern : the pattern of the split files in the disk preceding the incremental number, default="data_part"
#' @param type : a string c(rows, cols, list) indicating how to merge the data. This should be the same as the one 
#' used in the splitData procedure.
#' @param verbose : print running status information (default = TRUE).
#' @return a merged dataframe : Caution ! The merged results may be shuffeled. It is up to the user to reorder
#' the data accordingly.
mapResults <- function (folder = ".", pattern = "data_part", type = "rows", verbose=TRUE) 
{
  if (type %in% c("rows", "cols", "list") == FALSE)    # error corrected elc
  {
    stop("Please provide a valid mapping type rows, cols, list")
  }
  files <- dir(folder, pattern = paste(pattern, sep = "."))
  files.share <- as.numeric(gsub("data_part_|.rda", "", files))
  files <- files[order(files.share)]
  if (type %in% c("rows", "cols")) 
  {
    res <- readRDS(paste(folder, files[1], sep = "/"))
    for (i in 2:length(files)) 
    {
      if (i%%10 == 0 & verbose) 
      {
        print(paste("i =", i))
      }
      data_part_treated <- readRDS(paste(folder, files[i], sep = "/"))
      if (type == "rows") 
      {
        res <- rbind(res, data_part_treated) # error corrected elc
      } else 
      {
        res <- cbind(res, data_part_treated)
      }
      gc(verbose = FALSE)
    }
  } else 
  { # if LIST
    res <- list()
    for (i in 1:length(files)) 
    {
      if (i%%10 == 0 & verbose)  
      {
        print(paste("i =", i))
      }
      data_part_treated <- readRDS(paste(folder, files[i], sep = "/"))
      res <- c(res, data_part_treated) # error corrected elc
      gc(verbose = FALSE)
    }
  }
  warning("Caution! The merged results may be shuffeled. It is up to the user to reorder the data accordingly")
  return(res)
}


#' \code{deleteData} 
#' @title deleteData
#' @export
#' @description This function will delete the original and treated split temporary files to celan the workspace.
#' @author Edi Prifti
#' @param name : the name of the split files in the disk preceding the incremental number, default="data_part"
#' @return nothing to be returned
deleteData <- function(name="data_part"){
  files <- dir(pattern=paste(name,sep="_"))
  for(i in 1:length(files)){
    system(paste("rm -f",files[i]))
  }
}

