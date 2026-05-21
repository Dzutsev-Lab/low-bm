library(data.table)

#---------------------------
# Read nodes.dmp
#---------------------------
read_kraken_nodes <- function(nodes_file) {
  nodes_df <- fread(
    nodes_file,
    sep = "|",
    header = FALSE,
    fill = TRUE,
    quote = "",
    data.table = TRUE
  )

  nodes_df <- nodes_df[, .(
    TaxID = as.integer(trimws(V1)),
    Parent = as.integer(trimws(V2)),
    rank = trimws(V3)
  )]

  nodes_df <- nodes_df[!is.na(TaxID) & !is.na(Parent)]
  nodes_df
}

#---------------------------
# Build lookup tables
#---------------------------
build_parent_lookup <- function(nodes_df) {
  setNames(nodes_df$Parent, as.character(nodes_df$TaxID))
}

build_rank_lookup <- function(nodes_df) {
  setNames(nodes_df$rank, as.character(nodes_df$TaxID))
}

#---------------------------
# Taxonomy traversal helpers
#---------------------------
ancestor_chain <- function(taxid, parent_lookup) {
  taxid <- as.character(as.integer(taxid))
  if (is.na(taxid) || taxid == "0") return(integer(0))

  out <- integer(0)
  seen <- character(0)
  cur <- taxid

  while (!is.na(cur) && cur != "0") {
    if (cur %in% seen) break
    seen <- c(seen, cur)

    out <- c(out, as.integer(cur))

    next_parent <- parent_lookup[[cur]]
    if (is.null(next_parent) || is.na(next_parent)) break
    cur <- as.character(next_parent)
  }

  out
}

is_ancestor_of <- function(ancestor, descendant, parent_lookup) {
  ancestor <- as.integer(ancestor)
  descendant <- as.integer(descendant)

  if (is.na(ancestor) || is.na(descendant)) return(FALSE)
  if (ancestor == descendant) return(TRUE)

  ancestor %in% ancestor_chain(descendant, parent_lookup)
}

lowest_common_ancestor <- function(a, b, parent_lookup) {
  chain_a <- ancestor_chain(a, parent_lookup)
  chain_b <- ancestor_chain(b, parent_lookup)

  common <- intersect(chain_a, chain_b)
  if (length(common) == 0) return(0L)

  as.integer(common[1])
}


#-------------------------------------------
# Confidence score for assigned Kraken label
#-------------------------------------------
assigned_label_score <- function(kmer_trace, assigned_TaxID, parent_lookup) {
  if (is.na(kmer_trace) || is.na(assigned_TaxID)) return(NA_real_)

  tokens <- unlist(strsplit(trimws(kmer_trace), "\\s+"))
  if (length(tokens) == 0) return(NA_real_)

  C <- 0L
  Q <- 0L

  for (tok in tokens) {
    m <- regexec("^([^:]+):(\\d+)$", tok)
    parts <- regmatches(tok, m)[[1]]
    if (length(parts) != 3) next

    label <- parts[2]
    count <- as.integer(parts[3])
    if (is.na(count)) next

    if (label != "A") {
      Q <- Q + count
    }

    if (label != "A" && label != "0") {
      if (is_ancestor_of(assigned_TaxID, label, parent_lookup)) {
        C <- C + count
      }
    }
  }

  if (Q == 0L) return(NA_real_)
  C / Q
}