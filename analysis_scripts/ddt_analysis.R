pdf("ddt_analysis.pdf");

library("ggplot2");
library("reshape");

data <- read.table("ddtbench.out", header=TRUE);

## exclude some tests
#data <- data[data$testname != "FFT2",];
#data <- data[data$testname != "NAS_LU_x",];
##data <- data[data$testname != "NAS_LU_y",];
#data <- data[data$testname != "NAS_MG_x",];
#data <- data[data$testname != "NAS_MG_y",];
#data <- data[data$testname != "NAS_MG_z",];
data <- data[data$testname != "WRF_x_vec",];
data <- data[data$testname != "WRF_y_vec",];
##data <- data[data$testname != "LAMMPS_atomic",];
##data <- data[data$testname != "LAMMPS_full",];
##data <- data[data$testname != "SPECFEM3D_cm",];
data <- data[data$testname != "SPECFEM3D_oc",];
#data <- data[data$testname != "SPECFEM3d_mt",];
##data <- data[data$testname != "MILC_su3_zd",];

## helper functions

si_num <- function (x) {

    if (!is.na(x)) {
        if (x > 1e6) {
            chrs <- strsplit(format(x, scientific=12), split="")[[1]];
            rem <- chrs[seq(1,length(chrs)-6)];
            rem <- append(rem, "M");
        }

        else if (x > 1e3) {
            chrs <- strsplit(format(x, scientific=12), split="")[[1]];
            rem <- chrs[seq(1,length(chrs)-3)];
            rem <- append(rem, "K");
        }
        else {
        return(x);
    }
        retval <- paste(rem, sep="", collapse="");
        return(retval);
    }
    else return(NA);
}

si_vec <- function(x) {
  sapply(x, FUN=si_num);
}

ndata <- cast(data=data, formula=testname+bytes ~ packMethod+id, value="time", fun.aggregate=median);

## packing = pack+unpack
ndata$manual_packing = ndata$manual_pack + ndata$manual_unpack;

## build a new dataframe, first data for manual packing, then mpi ddts
testname <- as.character(ndata$testname);
bytes <- ndata$bytes;
value <- (ndata$bytes / ((ndata$manual_communication + ndata$manual_packing) / 2))
pp_data <- aggregate(data.frame(ndata$bytes, ndata$reference_communication),
                     by=list(Bytes = ndata$bytes), FUN=min)
testname <- append(testname, rep("Traditional Ping-Pong", length(pp_data[,1])));
bytes <- append(bytes, pp_data[,1]);
value <- append(value, (pp_data[,1] / (pp_data[,3] / 2)));
#packMethod <- rep("manual", length(bytes));
packMethod <- rep("mpi_pack_ddt", length(bytes));


testname <- append(testname, as.character(ndata$testname));
bytes <- append(bytes, ndata$bytes);
value <- append(value, ndata$bytes / (ndata$mpi_ddt_communication / 2));
pp_data <- aggregate(data.frame(ndata$bytes, ndata$reference_communication),
                     by=list(Bytes = ndata$bytes), FUN=min)
testname <- append(testname, rep("Traditional Ping-Pong", length(pp_data[,1])));
bytes <- append(bytes, pp_data[,1]);
value <- append(value, (pp_data[,1] / (pp_data[,3] / 2)));
#packMethod <- append(packMethod, rep("mpi_pack_ddt", length(ndata$bytes) + length(pp_data[,1])));
packMethod <- append(packMethod, rep("mpi_pack_ddt_dbrew", length(ndata$bytes) + length(pp_data[,1])));


newdf <- data.frame(testname, packMethod, bytes, value);

## From: http://stackoverflow.com/questions/38895918/value-missing-with-facet-labeller-function/38897092#38897092
labels_map <- c("manual" = "Manual Packing",
                "mpi_ddt" = "MPI DDTs",
                "mpi_pack_ddt" = "MPI DDT Pack/Unpack",
                "mpi_pack_ddt_dbrew" = "MPI DDT Pack/Unpack with DBrew",
                "mpi_pack_ddt_dbrew_llvm" = "MPI DDT P/U with DBrew/LLVM")

m <- ggplot(newdf, aes(x=bytes, y=value, color=testname, shape=testname)) +
    facet_grid(. ~ packMethod, scales="free", labeller=as_labeller(labels_map)) +
    geom_line() +
    geom_point(size=4) +
    scale_x_continuous('Datasize [Byte]', labels=si_vec) +
    scale_y_continuous('Bandwidth [MB/s]') +
    scale_colour_discrete(name="Test Name") +
    scale_shape_manual(name="Test Name", values=seq(1, length(unique(ndata$testname))+1)) +
    theme(axis.text.x = element_text(angle=90));

##print(m);
rm(testname, bytes, value, packMethod, pp_data, newdf, m)

## create new dataframe [testname  method  bytes  id   time] but
## for the manual method, the communication id is the sum of comm + pack + unpack

ndata <- cast(data=data, formula=testname+bytes+packMethod ~ id, value="time", fun.aggregate=median);

## we can do this for all methods, because for reference and mpi_ddt pack/unpack is 0.
ndata$communication <- ndata$communication + ndata$pack + ndata$unpack;

mndata <- melt(ndata, id=c("testname", "packMethod", "bytes"));

## do not plot reference and mpi_pack_ddt
mndata <- mndata[mndata$packMethod != "reference", ];
#mndata <- mndata[mndata$packMethod != "mpi_pack_ddt", ];
mndata <- mndata[mndata$packMethod != "mpi_ddt", ];
mndata <- mndata[mndata$packMethod != "manual", ];


## remove pack/unpack for mpi_ddt method, because it us 0 anyway
##mndata <- mndata[(mndata$packMethod != "mpi_ddt") | (mndata$id != "pack"),];
##mndata <- mndata[(mndata$packMethod != "mpi_ddt") | (mndata$id != "unpack"),];
                                        #
## Remove DBrew rewriting timing... Should be amortized ;)
mndata <- mndata[(mndata$packMethod != "mpi_pack_ddt") | (mndata$id != "dbrew_rewrite"),];
mndata <- mndata[(mndata$packMethod != "mpi_pack_ddt_dbrew") | (mndata$id != "dbrew_rewrite"),];
mndata <- mndata[(mndata$packMethod != "mpi_pack_ddt_dbrew_llvm") | (mndata$id != "dbrew_rewrite"),];

## We aren't interested in DDT create/free overheads at the moment.
mndata <- mndata[(mndata$id != "ddt_create_overhead"),];
mndata <- mndata[(mndata$id != "ddt_free_overhead"),];

print(mndata);


makeplot <- function(name) {
    print(name);
    m <- ggplot(mndata[mndata$testname == name,], aes(x=bytes, y=value, color=id, shape=id)) +
        geom_line() +
        facet_grid(. ~ packMethod, labeller=as_labeller(labels_map)) +
        geom_point(size=4) +
        labs(title=name) +
        scale_x_continuous('Datasize [Byte]', labels=si_vec) +
        scale_y_continuous('Time [us]') +
        scale_shape_manual(name="Phase", breaks=c("communication", "ddt_create_overhead", "ddt_free_overhead", "pack", "unpack", "dbrew_rewrite"),
                           labels=c("Communication RTT\n (incl. Packing, excl. DDT create)", "DDT Create", "DDT Free", "Pack", "Unpack", "Rewriting"),
                           values=c(21,22,23,24,25,26)) +
        scale_color_manual(name="Phase", breaks=c("communication", "ddt_create_overhead", "ddt_free_overhead", "pack", "unpack", "dbrew_rewrite"),
                           labels=c("Communication RTT\n (incl. Packing, excl. DDT create)", "DDT Create", "DDT Free", "Pack", "Unpack", "Rewriting"),
                           values=c(2:7)) +
        theme(axis.text.x = element_text(angle=90), legend.position="bottom");
    print(m);
    return(m);
}

plots = lapply(as.character(unique(mndata$testname)), makeplot);

dev.off()
