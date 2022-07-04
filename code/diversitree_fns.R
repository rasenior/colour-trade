pp.clades.fan <- diversitree:::pp.clades.fan
pp.clades.phylogram <- diversitree:::pp.clades.phylogram
pp.coords.fan <- diversitree:::pp.coords.fan
pp.coords.fix.xy <- diversitree:::pp.coords.fix.xy
pp.coords.phylogram <- diversitree:::pp.coords.phylogram
pp.lim.fan <- diversitree:::pp.lim.fan
pp.lim.phylogram <- diversitree:::pp.lim.phylogram
pp.node.coords <- diversitree:::pp.node.coords
pp.segments.fan <- diversitree:::pp.segments.fan
pp.segments.phylogram <- diversitree:::pp.segments.phylogram
pp.tiplabel.fan <- diversitree:::pp.tiplabel.fan
pp.tiplabel.phylogram <- diversitree:::pp.tiplabel.phylogram
group.label.tip <- diversitree:::group.label.tip
filled.arcs <- diversitree:::filled.arcs

plot2.phylo <-
    function(x, type = "phylogram", use.edge.length = TRUE, 
             node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE, 
             inner.branch.color = "black",
             edge.color = "black", edge.width = 1, edge.lty = 1, 
             font = 3, cex = par("cex"), adj = NULL, srt = 0, no.margin = FALSE, 
             root.edge = FALSE, label.offset = 0, underscore = FALSE, 
             x.lim = NULL, y.lim = NULL, direction = "rightwards", 
             lab4ut = "horizontal", tip.color = "black", ..., 
             n.taxa = NULL, clade.lwd = 1, clade.color = NULL, clade.fill = NULL, 
             pad = 0) {
        n.tip <- length(x$tip.label)
        if (n.tip < 2) stop("Cannot plot tree with < 2 tips")
        if (any(tabulate(x$edge[, 1]) == 1)) {
            stop("there are single (non-splitting) nodes in your tree;\n", 
                 "\tYou may need to use collapse.singles()")
        }
        type <- match.arg(type, c("phylogram", "fan"))
        direction <- match.arg(direction, "rightwards")
        if (is.null(x$edge.length)) use.edge.length <- FALSE
        if (no.margin) {
            par(mar = rep(0, 4))
        }
        tip.color <- rep(tip.color, length.out = n.tip)
        
        # edge.color <- rep(edge.color, length.out = nrow(x$edge))
        # Only colour branches that connect to a tip 
        # (i.e. where node is less than or equal to number of tips)
        edge.color <-
            sapply(1:nrow(x$edge), function(r){
                ifelse(x$edge[r, 2] <= length(x$tip.label), 
                       edge.color[x$edge[r, 2]], inner.branch.color)
            })
        edge.width <- rep(edge.width, length.out = nrow(x$edge))
        edge.lty <- rep(edge.lty, length.out = nrow(x$edge))
        if (!is.null(n.taxa)) {
            if (!is.null(names(n.taxa))) {
                n.taxa <- n.taxa[x$tip.label]
                n.taxa[is.na(n.taxa)] <- 1
            }
            else if (length(n.taxa) != n.tip) {
                stop("n.taxa must be same length as x$tip.label")
            }
        }
        x$n.taxa <- n.taxa
        x$n.spp <- length(x$tip.label) + sum(n.taxa - 1)
        if (!is.null(x$n.taxa)) {
            i <- match(seq_len(n.tip), x$edge[, 2])
            rep.clade <- function(x, edge) {
                if (is.null(x)) 
                    x <- edge
                else if (length(x) == 1) 
                    x <- rep(x, n.tip)
                else if (length(x) != n.tip) 
                    stop("Not yet handled")
                x
            }
            clade.color <- rep.clade(clade.color, edge.color[i])
            clade.fill <- rep.clade(clade.fill, edge.color[i])
            clade.lwd <- rep.clade(clade.lwd, edge.width[i])
            edge.color[match(which(x$n.taxa > 1), x$edge[, 2])] <- NA
        }
        pp.coords <- get(paste("pp.coords", type, sep = "."))
        pp.lim <- get(paste("pp.lim", type, sep = "."))
        pp.segments <- get(paste("pp.segments", type, sep = "."))
        pp.tiplabel <- get(paste("pp.tiplabel", type, sep = "."))
        pp.clades <- get(paste("pp.clades", type, sep = "."))
        asp <- ifelse(type == "fan", 1, NA)
        if (is.null(node.pos)) {
            node.pos <- 1
        } else if (node.pos != 1) {
            stop("node.pos != 1 not yet implemented")
        }
        if (is.null(adj)) adj <- 0
        xy <- pp.node.coords(x)
        xy.seg <- pp.coords(x, xy)
        lims <- pp.lim(x, xy, x.lim, y.lim, cex, show.tip.label, 
                       label.offset + pad)
        plot(NA, type = "n", xlim = lims$xlim, ylim = lims$ylim, 
             xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
             bty = "n", asp = asp,
             ...
             )
        pp.segments(x, xy.seg, edge.color, edge.width, edge.lty)
        if (!is.null(x$n.taxa) && any(x$n.taxa > 1)) {
            pp.clades(x, xy, xy.seg, clade.color, clade.fill, clade.lwd)
        }
        if (show.tip.label) {
            if (!underscore) 
                x$tip.label <- gsub("_", " ", x$tip.label)
            pp.tiplabel(x, xy, label.offset, adj, cex, tip.color, 
                        font)
        }
        if (type == "fan") {
            xy <- pp.coords.fix.xy(x, xy)
        }
        ret <- list(type = type, use.edge.length = use.edge.length, 
                    node.pos = node.pos, show.tip.label = show.tip.label, 
                    show.node.label = show.node.label, font = font, cex = cex, 
                    adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
                    x.lim = lims$xlim, y.lim = lims$ylim, direction = direction, 
                    tip.color = tip.color, Ntip = n.tip, Nnode = x$Nnode, 
                    edge = x$edge, xx = xy$xx, yy = xy$yy, n.taxa = x$n.taxa, 
                    n.spp = x$n.spp, xy = xy, xy.seg = xy.seg)
        assign("last_plot.phylo", ret, envir = .PlotPhyloEnv)
        invisible(ret)
    }


trait.plot <- 
    function(tree, dat, cols, lab = names(cols), str = NULL, class = NULL, 
              type = "f", w = 1/50, legend = length(cols) > 1, cex.lab = 0.5, 
              font.lab = 3, cex.legend = 0.75, margin = 1/4, check = TRUE, 
              quiet = FALSE, ...) 
    {
        if (!(type %in% c("f", "p"))) 
            stop("Only types 'f'an and 'p'hylogram are available")
        if (!is.null(class) && length(class) != length(tree$tip.label)) 
            stop("'class' must be a vector along tree$tip.label")
        n <- length(cols)
        if (n < 1) 
            stop("Need some colours")
        if (!is.data.frame(dat)) {
            if (is.vector(dat) && n == 1) {
                nm <- names(dat)
                dat <- matrix(dat)
                rownames(dat) <- nm
            }
            else {
                stop("dat must be a matrix")
            }
        }
        if (!all(tree$tip.label %in% rownames(dat))) 
            stop("All taxa must have entries in 'dat' (rownames)")
        if (n > 1) {
            if (!all(names(cols) %in% names(dat))) 
                stop("Not all colours have data")
            if (is.null(names(cols))) 
                stop("'cols' must be named")
            dat <- dat[names(cols)]
        }
        if (is.null(str)) {
            str <- lapply(dat, function(x) as.character(sort(unique(x))))
        }
        dat <- dat[tree$tip.label, , drop = FALSE]
        par(mar = rep(0, 4))
        t <- max(branching.times(tree))
        w <- w * t
        if (is.null(class)) {
            plt <- plot2.phylo(tree, type = type, show.tip.label = TRUE, 
                               label.offset = (n + 2) * w, cex = cex.lab, ...)
        }
        else {
            plt <- plot2.phylo(tree, type = type, show.tip.label = FALSE, 
                               label.offset = t * margin, ...)
            group.label.tip(plt, class, "black", "black", 
                            lwd = 1.5, offset.bar = w * (n + 2), offset.lab = w * 
                                (n + 3), cex = cex.lab, font = font.lab, check = check, 
                            quiet = quiet)
        }
        if (type == "f") {
            xy <- plt$xy
            theta <- xy$theta[seq_along(tree$tip.label)]
            dt <- diff(sort(theta))[1]/2
            for (i in seq_along(cols)) {
                idx <- dat[[names(dat)[i]]]
                if (any(idx == 0, na.rm = TRUE)) 
                    idx <- idx + 1
                filled.arcs(theta - dt, theta + dt, max(xy$x) + i * 
                                w, w, cols[[i]][idx])
            }
        }
        else {
            xy <- plt$xy[seq_along(tree$tip.label), ]
            dy <- 0.5
            for (i in seq_along(cols)) {
                idx <- dat[[names(dat)[i]]]
                if (any(idx == 0, na.rm = TRUE)) 
                    idx <- idx + 1
                xleft <- xy[1, 1] + w * i
                xright <- xleft + w
                ybottom <- xy[, 2] - dy
                ytop <- ybottom + dy * 2
                rect(xleft, ybottom, xright, ytop, col = cols[[i]][idx], 
                     border = NA)
            }
        }
        if (legend) {
            for (i in seq_along(cols)) {
                c.i <- cols[[i]]
                leg.txt <- str[[i]]
                leg.arg <- list(legend = leg.txt, title = lab[i], 
                                title.adj = 0, bty = "n", fill = c.i, cex = cex.legend, 
                                horiz = TRUE)
                ifelse(i == 1, leg <- do.call("legend", c("topleft", 
                                                          leg.arg)), leg <- do.call("legend", c(leg$rect$left, 
                                                                                                leg$rect$top - leg$rect$h, leg.arg)))
            }
        }
        invisible(plt)
    }