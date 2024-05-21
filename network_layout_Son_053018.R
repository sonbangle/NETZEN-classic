#!/usr/bin/Rscript
library(igraph)
library(optparse)
require(data.table)
require(plyr)
require(fastmatch)
#library(reticulate)
library(readr)
# install.packages("pheatmap", "RColorBrewer", "viridis")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(heatmaply)
library(optparse)
script_name = "network_layout_Son_053018.R"
source_folder = "$SOURCE"
data_folder = "${MY_BIN}/data"
seed = 0
set.seed(seed)

# get_full_source = function(source_file, source_folder=source_folder )
# {
#   source_path = system(paste("echo ", source_folder), intern = TRUE)
#   script = paste0(source_path,"/.., source_file)
#   return(script)
# }
source_path = system(paste("echo ", source_folder), intern = TRUE)
data_path = system(paste("echo ", data_folder), intern = TRUE)
#source_path = "/run/user/1000/gvfs/sftp:host=sftp.rc.ufl.edu/home/son.le/SOURCE_CODE/source"
#data_path = "/run/user/1000/gvfs/sftp:host=sftp.rc.ufl.edu//ufrc/dtran/son.le/MY_BIN/data"
#source_path = "/run/user/1835174/gvfs/sftp:host=sftp.rc.ufl.edu/home/son.le/SOURCE_CODE/source"
#data_path = "/run/user/1835174/gvfs/sftp:host=sftp.rc.ufl.edu//ufrc/dtran/son.le/MY_BIN/data"


MSigDB_file = paste0(data_path, "/msigdb_v6.1.xml")
source(paste0(source_path, "/utils.R"))
source(paste0(source_path, "/display_heatmap.R"))
job_log_dir = "./cluster_logs/jobs" # Folder containing job files to submit to computer cluster. Global variable
job_count = 0 # Index of job. Use in display_heatmap.R or any place that need to submit job to cluster. job script will be placed in job_log_dir folder. Global variable


excluded_word_list = c(
  "from",
  "to",
  "of",
  "regulation",
  "binding",
  "activity",
  "activated",
  "protein",
  "complex",
  "positive",
  "negative",
  "process",
  "signaling",
  "response",
  "pathway",
  "involved",
  "by",
  "in",
  "into",
  "across",
  "NA",
  "production",
  "biosynthetic",
  "factor" ,
  "gene",
  "genes",
  "Genes",
  "the",
  "upregulated",
  "downregulated",
  "and",
  "up-regulated",
  "down-regulated",
  "Neighborhood",
  "comparison",
  "with",
  "versus",
  "after",
  "at",
  "module",
  "that",
  "compared",
  "any",
  "around",
  "a",
  "control",
  "disease" ,
  "after",
  "before",
  "cells",
  "top",
  "effect",
  "expression",
  "gender",
  "global",
  "profiling",
  "-",
  "activation",
  "through",
  "contain",
  "off",
  "is",
  "This",
  "this",
  "The",
  "RefSeq,",
  "[provided",
  "encoded",
  "encodes",
  "for",
  "which",
  "been",
  "an",
  "be",
  "have",
  "are",
  "family",
  "may",
  "on",
  "but",
  "role",
  "encoding",
  "as",
  "2008]",
  "June",
  "Jul",
  "Dec",
  "has",
  "transcript",
  "variants",
  "proteins",
  "different",
  "member",
  "found",
  "results",
  "composed",
  "functions",
  "also",
  "It",
  "Core",
  "gene.",
  "isoforms",
  "approximately",
  "family.",
  "play",
  "function",
  "other",
  "identified",
  "(MIM",
  "one",
  "isoforms",
  "expressed",
  "multiple",
  "splicing",
  "Alternative",
  "including",
  "human",
  "some",
  "cell",
  "associated",
  "spliced",
  "type",
  "Mutations",
  "its",
  "several",
  "subunit",
  "subunit.",
  "variants.",
  "Alternatively",
  "contains",
  "domain",
  "enzyme",
  "or",
  "form",
  "belongs",
  "receptor",
  "cluster",
  "members",
  "receptors",
  "two",
  "belongs",
  "not",
  "containing",
  "like",
  "1",
  "2",
  "3",
  "4",
  "5",
  "6",
  "7",
  "8",
  "9",
  "0",
  "molecule",
  "motif",
  "band",
  "cytogenetic"
  
  
  
  
  
)

excluded_GO_terms = c(
  "signal transduction",
  "molecular_function",
  "biological_process",
  "cellular_component",
  "membrane",
  "extracellular region",
  "integral component of plasma membrane",
  "extracellular space",
  "extracellular matrix",
  "intracellular",
  "cell surface",
  "protein complex",
  "protein transport",
  "metal ion binding",
  "regulation of gene expression",
  "sequence-specific DNA binding",
  "transcription factor binding",
  "protein homodimerization activity",
  "positive regulation of transcription from RNA polymerase II promoter "
  ,
  "positive regulation of transcription from RNA polymerase II promoter",
  "protein homodimerization activity",
  "protein heterodimerization activity" ,
  "Systemic lupus erythematosus"
)




get_communities = function(graph, method)
{
  if (method == "louvain")
  {
    communities = cluster_louvain(graph, weight = NA)
    
  } else
    if (method == "infomap")
    {
      communities = cluster_infomap(graph)
    }
  return(communities)
}


get_subnets = function(graph, clustering_method)
{
  communities = get_communities(graph, clustering_method)
  n_communities = length(communities)
  if (n_communities == 1)
  {
    #   print("28, n  communities =1, now use louvain to breakdown graph")
    #   communities = get_communities(graph,method="louvain")  # if infomap cannot breakdown a graph further, then switch to louvain as louvain can breakdown any grapy
    #
    #   n_communities = length(communities)
    print(" only one community, cannot breakdown further")
  }
  
  n_nodes = communities$vcount
  communities_list  = list()
  membership = communities$membership
  for (i in c(1:n_communities))
  {
    communities_list[[i]] = list()
  }
  for (i in c(1:n_nodes))
  {
    cluster = membership[i]
    communities_list[[cluster]] = c(communities_list[[cluster]], i)
  }
  n_communities = length(communities)
  subnets = list()
  for (i in c(1:n_communities))
  {
    community = communities_list[[i]]
    subnet = induced_subgraph(graph, community)
    subnets[[i]] = subnet
  }
  return(list(subnets, n_communities))
}

# attributes_list
get_hierarchy_clusters = function(subnets,
                                  clustering_method = "louvain",
                                  max_cluster_size = 100,
                                  ndim = 3,
                                  layout_method = "fr",
                                  level = 1,
                                  prefix = "",
                                  subnets_layout,
                                  space_ratio = 2,
                                  node_score_field = NULL,
                                  beta = 0.5,
                                  parent_score = 0,
                                  size_field = "degree",
                                  attributes_list = None,
                                  parent_subnet_attributes_list = None)
{
  # beta is the coefficient for determining the node final score when doing cluster smoothing.
  # attributes is a data frame with columns(gene, attributues) where SYMBOL is gene symbol, attributes is attribute names (such as Gen Ontology ID), rownames is gene symbol
  subnet_list = list()
  if (!is.null(attributes_list) &
      class(attributes_list) != "list")
  {
    attributes_list = list(attributes_list)
  }
  cluster_radius = get_radius(
    subnets = subnets,
    subnets_layout = subnets_layout,
    space_ratio = space_ratio,
    ndim = ndim,
    size_field = size_field,
    net_type = 1
  )
  col_axis = c()
  if (!is.null(attributes_list))
  {
    for (i in c(1:length(attributes_list)))
    {
      attributes = attributes_list[[i]]
      rownames(attributes) = attributes[, 1]
      attributes_list[[i]] = attributes
    }
    
  }
  for (i in c(1:ndim))
  {
    col_axis = c(col_axis, paste("X", i, sep = ""))
  }
  for (i in c(1:length(subnets)))
  {
    subnet = subnets[[i]]
    n_node = length(V(subnet))
    element_node_names = V(subnet)$name
    #print(paste(i, "59", n_node, max_cluster_size))
    center = subnets_layout[i, ]
    radius = cluster_radius[[i]]
    if (prefix != "")
    {
      subnet_name = paste(prefix, ".", i, sep = "")
    }
    else
    {
      subnet_name = as.character(i)
    }
    print(paste("63", "subnet name:", subnet_name, "radius:", radius))
    if (radius > 1)
    {
      print(paste(82, radius))
      break()
      
    }
    if (!is.null(node_score_field))
    {
      subnet_score = vertex_attr(subnet, name = node_score_field)
      subnet_score[is.na(subnet_score)] = 0
      subnet_score_genes = subnet_score
      subnet_score = mean(abs(subnet_score))
      subnet_score = beta * subnet_score + (1 - beta) * parent_score  # pull the subnet score to the parent score
      
      #Calculate change direction (positive or negative, overexpressed or downregulated) by number of genes that over express - number of genes that down regulated, if that >0: then overexpress)
      n_up = length(subnet_score_genes[subnet_score_genes > 0])
      n_down = length(subnet_score_genes[subnet_score_genes < 0])
      
      #change = mean(subnet_score_genes)
      change = median(subnet_score_genes[subnet_score_genes != 0])
      
      print(paste("378 network_layout_Son, change", change))
      if (is.na(change))
      {
        print("change is NULL")
        print(subnet_score_genes)
        print(158)
        print(change)
        print(subnet_score)
        print(subnet)
        print(node_score_field)
        print(vertex_attr(subnet))
        change = 0
      }
      if (change > 0)
      {
        subnet_change_direction = 1
      }
      if (change == 0)
      {
        subnet_change_direction = 0
      }
      if (change < 0)
      {
        subnet_change_direction = -1
      }
      
      
    } else
    {
      subnet_score = 0
      subnet_change_direction = 0
    }
    
    if (is.null(subnet_score))
    {
      break()
    }
    
    if (!is.null(attributes))
    {
      subnet_attributes_list = list()
      for (i in c(1:length(attributes_list)))
      {
        attributes = attributes_list[[i]]
        subnet_attributes  = attributes[fmatch(element_node_names, rownames(attributes)), ]
        #my.data[match(random.rows,rownames(my.data)),]
        subnet_attributes = colMeans(subnet_attributes[, c(2:ncol(subnet_attributes))], na.rm = TRUE)
        subnet_attributes[is.na(subnet_attributes)] = 0
        subnet_attributes = beta * subnet_attributes + (1  - beta) * parent_subnet_attributes_list[[i]]  # pull the subnet score to the parent score
        subnet_attributes_list[[i]] = subnet_attributes
      }
    }
    else
    {
      subnet_attributes_list = NULL
    }
    
    get_gene_layout = FALSE
    
    if (n_node > max_cluster_size)
    {
      old_subnet = subnet
      element_node_names = V(subnet)$name
      get_subnets_result = get_subnets(subnet, clustering_method)
      subnet = get_subnets_result[[1]]
      n_communities = get_subnets_result[[2]]
      if (n_communities > 1)
      {
        cluster_net = get_cluster_net(old_subnet, subnet, subnet_name)
        
        net_type = 1   # intermediate subnet containing sub subnets, not final subnet that contains only genes
        subnet_list_layout = update_subnet_list(
          subnet_list,
          cluster_net,
          subnet_name,
          layout_method,
          ndim,
          center,
          radius,
          subnet_score,
          n_genes = n_node,
          level = level,
          size_field = size_field,
          net_type = net_type,
          subnet = subnet,
          space_ratio = space_ratio,
          element_node_names = element_node_names,
          n_up = n_up,
          n_down = n_down,
          subnet_change_direction
        )
        subnet_list = subnet_list_layout[[1]]
        cluster_layout = subnet_list_layout[[2]]
        #print(paste("81 continue breakdown", n_node, max_cluster_size))
        #print(cluster_layout)
        subnet_result = get_hierarchy_clusters(
          subnet,
          clustering_method ,
          max_cluster_size,
          ndim = ndim,
          layout_method = layout_method,
          level = level + 1,
          prefix = subnet_name,
          subnets_layout = cluster_layout[, col_axis],
          space_ratio = space_ratio,
          node_score_field = node_score_field,
          beta = beta,
          parent_score = subnet_score,
          attributes_list = attributes_list,
          parent_subnet_attributes_list = subnet_attributes_list
        )
        for (net in (subnet_result[[1]]))
        {
          subnet_list[[length(subnet_list) + 1]] = net
        }
        layout = subnet_result[[2]]
        if (!is.null(attributes_list))
        {
          subnet_attributes_list = subnet_result[[3]]
        }
      }
      else
      {
        subnet = old_subnet
        get_gene_layout = TRUE
      }
      
    } else
    {
      # get layout
      get_gene_layout = TRUE
    }
    
    if (get_gene_layout)
    {
      net_type = 2
      subnet_list_layout = update_subnet_list(
        subnet_list,
        subnet,
        subnet_name,
        layout_method,
        ndim,
        center,
        radius,
        subnet_score,
        n_genes = n_node,
        level = level,
        size_field = size_field,
        net_type = net_type,
        subnet = subnet,
        space_ratio = space_ratio,
        element_node_names = element_node_names,
        n_up = n_up,
        n_down = n_down,
        subnet_change_direction
      )
      subnet_list = subnet_list_layout[[1]]
      layout = subnet_list_layout[[2]]
      layout = data.frame(layout, check.names = FALSE)
      layout$subnet_name = subnet_name
      layout$level = level + 1
      layout$net_type = 3
      # cluster_radius = get_radius(subnets = subnets, subnets_layout = subnets_layout, space_ratio = space_ratio, ndim = ndim)
      if (!is.null(node_score_field))
      {
        v_attr  = vertex_attr(subnet, name = node_score_field)
        v_attr[is.na(v_attr)] = 0
        v_attr = beta * abs(v_attr) + (1 - beta) * subnet_score
        layout$node_smooth_score = v_attr
        layout$subnet_score = subnet_score
      }
      if (!is.null(attributes_list))
      {
        for (i in c(1:length(attributes_list)))
        {
          attributes = attributes_list[[i]]
          subnet_attributes = subnet_attributes_list[[i]]
          attr_names = colnames(attributes)[2:ncol(attributes)]
          subnet_attributes = data.frame(t(subnet_attributes), check.names = FALSE)
          subnet_attributes = subnet_attributes[rep(1, nrow(layout)), ]
          new_subnet_attributes = data.frame(SYMBOL = element_node_names)
          subnet_attributes = cbind(new_subnet_attributes, subnet_attributes)
          subnet_attributes_list[[i]] = subnet_attributes
        }
        
      }
    }
    
    
    
    if (exists("end_layout"))
    {
      # if (ncol(layout) != ncol(end_layout))
      # {
      #   print("wrong")
      #   break()
      # }
      end_layout = rbind(end_layout, layout)
    } else
    {
      end_layout = layout
    }
    
    if (exists("end_attributes_list"))
    {
      # if (ncol(end_attributes) != ncol(subnet_attributes))
      # {
      #   print("wrong")
      #   break()
      # }
      
      for (i in c(1:length(attributes_list)))
      {
        end_attributes_list[[i]] = rbind(end_attributes_list[[i]], subnet_attributes_list[[i]])
        
        
        
      }
      
      
      
    } else
    {
      end_attributes_list = list()
      for (i in c(1:length(attributes_list)))
      {
        end_attributes_list[[i]] = subnet_attributes_list[[i]]
      }
      
    }
    
  }
  if (!is.null(attributes_list))
  {
    out = list(subnet_list, end_layout, end_attributes_list)
  } else
  {
    out = list(subnet_list, end_layout)
  }
  return(out)
}

# concatenate cluster_n
# n_up: number of upregulated genes
# n_down: number of downregulated genes
update_subnet_list = function(subnet_list,
                              cluster_net,
                              cluster_net_name,
                              layout_method,
                              ndim,
                              center,
                              radius,
                              subnet_score,
                              n_genes,
                              level,
                              size_field,
                              net_type,
                              subnet,
                              space_ratio,
                              element_node_names,
                              n_up = 0,
                              n_down = 0,
                              subnet_change_direction)
{
  cluster_layout = get_layout(cluster_net, layout_method, ndim, center, radius)
  
  #cluster_net = set_graph_attr(cluster_net, name="radius", value = radius)
  cluster_net = set_graph_attr(cluster_net, name = "graph_name", value = cluster_net_name)
  cluster_net = set_graph_attr(cluster_net, name = "n_genes", value = n_genes)
  cluster_net = set_graph_attr(cluster_net, name = "level", value = level)
  cluster_net = set_graph_attr(cluster_net, name = "net_type", value = net_type)
  cluster_net = set_graph_attr(cluster_net, name = "element_node_names", value = element_node_names)
  # subnet_change_direction = 0
  # try({
  #   if (n_up >= n_down)
  #   {
  #     subnet_change_direction= 1
  #   } else
  #   {
  #     subnet_change_direction= -1
  #   }
  #
  # })
  cluster_net = set_graph_attr(cluster_net, name = "change_direction", value = subnet_change_direction)
  cluster_net = set_graph_attr(cluster_net, name = "n_up", value = n_up)
  cluster_net = set_graph_attr(cluster_net, name = "n_down", value = n_down)
  for (i in c(1:ndim))
  {
    axis = paste("X", i, sep = "")
    cluster_net = set_graph_attr(cluster_net, name = axis, value = center[1, axis])
    cluster_net = set_vertex_attr(cluster_net, name = axis, value = cluster_layout[, axis])
  }
  if (!is.null(subnet_score))
  {
    cluster_net = set_graph_attr(cluster_net, name = "subnet_score", value = subnet_score)
  }
  cluster_radius = get_radius(
    subnets = subnet,
    subnets_layout = cluster_layout,
    space_ratio = space_ratio,
    ndim = ndim,
    size_field = size_field,
    net_type = net_type
  )
  cluster_radius = unlist(cluster_radius)
  cluster_net = set_vertex_attr(cluster_net, name = "radius", value = cluster_radius)
  cluster_layout$radius = cluster_radius
  subnet_list[[length(subnet_list) + 1]] = cluster_net
  # for (i in c(1:n(cluster_layout)))
  for (i in c(1:nrow(cluster_layout)))
  {
    for (j in c(1:ncol(cluster_layout)))
    {
      if (is.na(cluster_layout[i, j]))
      {
        print("NA")
      }
      
    }
  }
  
  return(list(subnet_list, cluster_layout))
}





# get_hierarchy_clusters_from_graph = function(graph, clustering_method, max_cluster_size, ndim, layout_method,space_ratio, node_score_field, beta)
# {
#   print("138 getting hierarchy clusters from graph ...")
#   print("139 getting first level subnets. ...")
#   subnets = get_subnets(graph, clustering_method)
#   print("first level subnets received")
#   merged_cluster_net = get_cluster_net (graph, subnets, prefix="1")
#   layout = get_layout(merged_cluster_net,method = "fr", ndim, center=, radius =1)
#   final_subnets = get_hierarchy_clusters(subnets, clustering_method,max_cluster_size,ndim = ndim, layout_method = layout_method,
#                                          subnets_layout = layout[,2:ncol(layout)], space_ratio = space_ratio,node_score_field = node_score_field, beta = beta, parent_score = 0.)
#   return(final_subnets)
# }



get_hierarchy_clusters_from_graph = function(graph,
                                             clustering_method = "louvain",
                                             max_cluster_size = 100,
                                             ndim = 3,
                                             layout_method = "fr",
                                             space_ratio,
                                             node_score_field,
                                             beta = 0.5,
                                             attr = None)
{
  subnets = list(graph)
  subnets_layout = data.frame()
  for (i in c(1:ndim))
  {
    subnets_layout[1, paste("X", i, sep = "")] = 0
  }
  
  if (!is.null(attr))
  {
    if (class(attr) != "list")
    {
      attr = list(attr)
    }
    parent_subnet_attributes_list = list()
    for (i in c(1:length(attr)))
    {
      att = attr[[i]]
      parent_subnet_attributes_list[[i]] = rep(0, ncol(att) - 1)
    }
    
  } else
  {
    parent_subnet_attributes_list = NULL
  }
  
  final_subnets = get_hierarchy_clusters(
    subnets,
    clustering_method,
    max_cluster_size,
    ndim = ndim,
    layout_method = layout_method,
    subnets_layout = subnets_layout,
    space_ratio = space_ratio,
    node_score_field = node_score_field,
    beta = beta,
    parent_score = 0.,
    prefix = "",
    attributes_list = attr,
    parent_subnet_attributes_list = parent_subnet_attributes_list
  )
  return(final_subnets)
}



get_layout = function(graph,
                      method = "fr",
                      ndim,
                      center = rep(0, ndim),
                      radius = 1)
{
  if (length(V(graph)) == 2)
    # graph have only two nodes
  {
    layout = data.frame(X1 = c(0, 1), X2 = c(0, 1))
    if (ndim == 3)
    {
      layout$X3 = c(0, 1)
    }
  } else
  {
    if (method == "fr")
    {
      layout = layout_with_fr(graph, dim =  ndim)
    }
    if (method == "kk")
    {
      layout = layout_with_kk(graph, dim =  ndim)
    }
    
    layout = data.frame(layout)
  }
  layout = normalize_layout(layout, center =  center, radius = radius)
  node_name = as.character(V(graph)$name)
  nodes = data.frame(node_name = node_name)
  layout = cbind(nodes, layout)
  for (i in c(1:nrow(layout)))
  {
    for (j in c(1:ncol(layout)))
    {
      if (is.na(layout[i, j]))
      {
        print("NA")
      }
      
    }
  }
  return(layout)
}


# get cluster net by merging nodes within subclusters, and create edges between subclusters
# input: net: a whole net containing all nodes and all edges
# subnets: a subnet of each cluster
get_cluster_net = function(net,
                           subnets,
                           prefix = "1",
                           edge_weight_field = NA)
{
  merged_net = make_empty_graph(directed = FALSE)
  
  n_cluster = length(subnets)
  merged_net = add_vertices(merged_net, n_cluster)
  if (prefix != "")
  {
    V(merged_net)$name = paste(prefix, ".", c(1:n_cluster), sep = "")
  } else
  {
    V(merged_net)$name = as.character(c(1:n_cluster))
  }
  node_subnet_dict = new.env(hash = TRUE, size = length(E(net)))
  
  for (cluster_index in c(1:n_cluster))
  {
    subnet_nodes = V(subnets[[cluster_index]])$name
    # print(paste(cluster_index, length(subnet_nodes)))
    for (node in subnet_nodes)
    {
      node_subnet_dict[[node]] = cluster_index
    }
    
    
  }
  
  
  if (is.na(edge_weight_field))
  {
    for (i in c(1:length(V(net))))
    {
      nodes = ends(net, i)
      source_cluster = node_subnet_dict[[nodes[1]]]
      target_cluster = node_subnet_dict[[nodes[2]]]
      if (source_cluster != target_cluster)
      {
        if (are_adjacent(merged_net, source_cluster, target_cluster))
        {
          edge_id = get.edge.ids(merged_net, c(source_cluster, target_cluster))
          E(merged_net)[[edge_id]]$weight = E(merged_net)[[edge_id]]$weight +
            1 / 100
        } else
        {
          # print(paste(171, source_cluster, target_cluster))
          merged_net = merged_net + edge(source_cluster, target_cluster, weight = 1 /
                                           100)
          #print(E(merged_net))
        }
      }
    }
  }
  #print(merged_net)
  
  return(merged_net)
}




# normalize layout so that the center of layout is the center, and the max distance from any point to the center is not more than radius * oversize_factor
normalize_layout = function(layout,
                            center = c(1., 1.),
                            radius = 1.,
                            oversize_factor = 1)
{
  ncols = ncol(layout)
  old_center = sapply(layout, mean)
  dist_to_center = apply(layout, 1, function(x)
  {
    dist = x - old_center
    dist = dist ^ 2
    dist = mean(dist)
    dist = sqrt(dist)
    return(dist)
  })
  max_dist_to_center = max(dist_to_center) * oversize_factor
  scale = radius / max_dist_to_center
  layout = layout * scale
  old_center = sapply(layout, mean)
  center = data.frame(center)
  translation = t(center - old_center)
  new_layout = apply(layout, 1, function(x)
  {
    return (x + translation)
  })
  new_layout = data.frame(t(new_layout))
  
  return(new_layout)
  
}

# Get appropirate radius for each cluster sphere so that the clusters do not overlap each other
# space ratio is to regulate the space between clusters, the more the space ratio, the more gap between clusters
# subnets can be a list of subnet graphs or one final subnet with nodes are individual genes)
get_radius = function(subnets,
                      subnets_layout,
                      space_ratio,
                      ndim,
                      size_field,
                      net_type = 2)
{
  if (net_type != 2)
  {
    n_subnets = length(subnets)
  }
  else
  {
    n_subnets = length(V(subnets))
  }
  cluster_radius = list()
  if (n_subnets == 1)
  {
    cluster_radius = list(1.)
    return(cluster_radius)
  }
  
  dist_ratios = list()
  if (net_type != 2)
    # subnets is not a list of gene node
  {
    for (i in c(1:n_subnets))
    {
      subnet = subnets[[i]]
      
      # cluster_radius[[i]] = length(V(subnet)) # simple as the degree of node
      
      size = vertex_attr(subnet, name = size_field)
      radius = size ^ (1 / ndim)
      radius  = sum(radius)
      
      cluster_radius[[i]] = radius
      
    }
  }
  else
  {
    radius_list = vertex_attr(subnets, name = size_field)
    for (i in c(1:n_subnets))
    {
      radius = radius_list[i]
      radius = radius ^ (1 / ndim)
      cluster_radius[[i]] = radius
    }
    
  }
  
  #print("354")
  #print(cluster_radius)
  
  axis_cols = colnames(subnets_layout)[2:ncol(subnets_layout)]
  col_axis = c()
  for (i in c(1:ndim))
  {
    col_axis = c(col_axis, paste("X", i, sep = ""))
  }
  
  nodes_layout = subnets_layout[, col_axis]
  dist_object = dist(nodes_layout)
  dist_matrix = as.matrix(dist_object)
  
  
  for (i in c(1:n_subnets))
  {
    for (j in c(1:n_subnets))
    {
      if (i != j)
      {
        #         pos_i = nodes_layout[i,]
        #         pos_j = nodes_layout[j,]
        #         dist = (pos_j - pos_i)
        #         dist = dist^2
        #         dist = sum(dist)
        #         dist = sqrt(dist)
        dist = dist_matrix[i, j]
        
        #print(paste(395, dist, dist2, dist2-dist, (dist2- dist)/dist))
        min_dist = cluster_radius[[i]] + cluster_radius[[j]]
        dist_ratio = min_dist / dist  # should be less than 1 to avoid collision. dist is fixed
        dist_ratios = c(dist_ratios, dist_ratio)
        #print(paste(403, dist, min_dist, dist_ratio))
      }
    }
  }
  max_dist_ratio = max(unlist(dist_ratios))
  scale = 1 / (max_dist_ratio * space_ratio)
  for (i in c(1:n_subnets))
  {
    # print(paste(411,  cluster_radius[[i]], scale ))
    cluster_radius[[i]] = cluster_radius[[i]] * scale
    # print(paste(413,  cluster_radius[[i]]))
    if (is.na(cluster_radius[[i]]))
    {
      print(387)
      break()
    }
    if (cluster_radius[[i]] > 10)
    {
      print(387)
      break()
    }
  }
  
  
  for (radius in cluster_radius)
  {
    if (is.na(radius) | radius > 1)
    {
      print(radius)
      break()
    }
  }
  
  return(cluster_radius)
}



get_network_layout = function(graph,
                              clustering_method,
                              max_cluster_size,
                              ndim,
                              layout_method,
                              space_ratio,
                              node_score_field,
                              beta,
                              dist_normalize = FALSE,
                              attr = NULL)
{
  print("295 getting network layout ..")
  clusters = get_hierarchy_clusters_from_graph(
    graph,
    clustering_method,
    max_cluster_size,
    ndim = ndim,
    layout_method = layout_method,
    space_ratio,
    node_score_field = node_score_field,
    beta = beta,
    attr = attr
  )
  final_subnet = clusters[[1]]
  final_layout = clusters[[2]]
  if (!is.null(attr))
  {
    node_attr_score = clusters[[3]]
  }
  if (dist_normalize == TRUE)
  {
    # now normalize the final layout so that the center of layout is [0, 0, 0] and max distance to the center  is 1.
    dims = paste("X", c(1:ndim), sep = "")
    locs = final_layout[, dims]
    
    center = sapply(locs, mean)
    #print(paste(318, "center", center))
    dists = apply(locs, 1, function(x)
    {
      dist = 0
      for (i in c(1:ndim))
      {
        axis_dist = x[i] - center[i]
        axis_dist = axis_dist ^ 2
        dist = dist + axis_dist
      }
      dist = sqrt(dist)
    })
    
    max_dist = max(dists)
    scale = 1 / max_dist
    final_layout[, dims] = final_layout[, dims] * scale
    new_center = center * scale
    for (i in c(1:ndim))
    {
      final_layout[, paste("X", i, sep = "")] = final_layout[, paste("X", i, sep =
                                                                       "")] - new_center[i]
    }
    locs = final_layout[, dims]
    
    center = sapply(locs, mean)
    #print(paste(343, "center", center))
    
    node_degree = unlist(degree(graph))
    node_names = V(graph)$name
    node_info = data.frame(name = node_names, degree = node_degree)
    final_layout = data.table(final_layout, key = "node_name")
    final_layout$row_ID = 1:nrow(final_layout)
    node_info = data.table(node_info, key = "name")
    final_layout = fast_merge(
      final_layout,
      node_info,
      by = c("node_name", "name"),
      all.x = TRUE,
      sort = FALSE
    )
    final_layout = final_layout[order(final_layout$id),]
    #print(426)
    #print(final_layout[1,])
    new_final_subnet = list()
    for (i in c(1:length(final_subnet)))
    {
      cluster = final_subnet[[i]]
      radius = graph_attr(cluster, "radius") * scale
      cluster = set_graph_attr(cluster, "radius", radius)
      for (i in c(1:ndim))
      {
        axis = paste("X", i, sep = "")
        old_coord = vertex_attr(cluster, name = axis)
        new_coord = old_coord * scale - new_center[i]
        cluster = set_vertex_attr(cluster, name = axis , value =  new_coord)
        old_subnet_center = graph_attr(cluster, axis)
        new_subnet_center = old_subnet_center * scale  - unlist(new_center[i])
        cluster = set_graph_attr(cluster, axis, new_subnet_center)
        graph_name = graph_attr(cluster)$graph_name
        
      }
      
      
      
      new_final_subnet[[length(new_final_subnet) + 1]] = cluster
    }
    final_subnet = new_final_subnet
    
  }
  
  if (is.null(attr))
  {
    return(list(final_layout, final_subnet))
  }
  else
  {
    {
      return(list(final_layout, final_subnet, node_attr_score))
    }
  }
  
  
}


get_enrichment_score = function(inputs, sampling_size = 1000)
{
  inputs = log(inputs + 0.0000001)
  random_data = sample(inputs, size = sampling_size)
  random_data = random_data[random_data > -6]
  random_mean = mean(random_data)
  random_sd = sd(random_data)
  p_value = pnorm(inputs,
                  mean = random_mean,
                  sd = random_sd,
                  lower.tail = FALSE)
  fdr = p.adjust(p_value, method = "BH")
  enrichment_score = (1 - fdr)
  return(enrichment_score)
}

# Similar to get_radius, get appropirate node_size for each node in graph so that nodes do not overlap each other
# space ratio is to regulate the space between nodes, the more the space ratio, the more gap between clusters
#layout should contain at least the following colums: node_name, degree, X1, X2, ...
# Return the same dataframa from layout but with additional column: node_radius
# size_field: field where to get value to calculate node radius, defaut is degree
# add_node_radius = function(layout, space_ratio, ndim, size_field = "degree")
# {
#   n_nodes = nrow(layout)
#   nodes_radius = layout[, size_field]^(1/ndim)
#
#   for (k in c(1:ndim))
#   {
#     if (k == 1) {axis = "X1"} else
#     {axis = c(axis, paste("X", k, sep=""))}
#   }
#
#   for (i in c(1:n_nodes))
#   {
#     print(paste(i, "calculating node size for ", layout[i,"node_name"]))
#     for (j in c(1:n_nodes))
#     {
#       if(i != j)
#       {
#         pos_i = layout[i,axis]
#         pos_j = layout[j,axis]
#         dist = (pos_j - pos_i)
#         dist = dist^2
#         dist = mean(dist)
#         dist = sqrt(dist)
#         min_dist = nodes_radius[i] + nodes_radius[j]
#         dist_ratio = dist / min_dist
#         if (exists("min_dist_ratio"))
#         {
#           if (dist_ratio < min_dist_ratio) { min_dist_ratio = dist_ratio}
#         } else min_dist_ratio = dist_ratio
#
#       }
#     }
#   }
#
#   scale = space_ratio / min_dist_ratio
#   layout$node_radius = nodes_radius /scale
#
#
#   return(layout)
# }




# The main pipeline from the input network file to the output for importing into vtk script
# node_attributes is a data frame with first column containing gene name and following columns contain attribute values for each gene
get_network_layout_pipeline = function(network_file = "../data/nha_gsc_net_nbhscore.graphml" ,
                                       node_attributes = NULL,
                                       graph_format = "graphml",
                                       #ncol
                                       clustering_method = "louvain",
                                       max_cluster_size = 100 ,
                                       ndim = 2,
                                       layout_method = "fr",
                                       space_ratio = 2,
                                       node_score_field = "logFC",
                                       enrichment_score_field = c("node_smooth_score"),
                                       sampling_size = 1000,
                                       beta = 0.5,
                                       n_top_annotation = 3,
                                       importance_field = "score",
                                       # Use for annotation. Problem is nscore result is not available for all nodes, and sometimes consider positive only node.
                                       n_top_genes = 3,
                                       outdir = "../data/layout",
                                       debug = TRUE)
{
  dir.create(outdir)
  if (debug == TRUE)
    
  {
    print("saving debug RDS Rdata files")
    inputs = list(
      network_file = network_file,
      node_attributes = node_attributes,
      graph_format = graph_format,
      clustering_method = clustering_method,
      max_cluster_size = max_cluster_size
    )
    debug_dir = paste0(outdir, "/layout/", space_ratio)
    dir.create(debug_dir)
    saveRDS(inputs,
            file = paste0(debug_dir, "/get_network_layout_pipeline.RDS"))
    save.image(file = paste0(debug_dir, "/get_network_layout_pipeline.Rdata"))
  }
  print("reading graph file ...")
  graph = read_graph(file = network_file, format = graph_format)
  
  print("graph file has been read")
  graph = as.undirected(graph, mode = c("collapse"))
  graph_degree = degree(graph)
  graph = set_vertex_attr(graph, name = "degree", value = graph_degree)
  nodes = V(graph)$name
  nodes = toupper(nodes) # Change node name to upper case to make sure consistency accross all networks
  V(graph)$name = nodes
  nodes = data.frame(SYMBOL = nodes)
  if (!is.null(node_attributes))
  {
    node_attributes = import_file(node_attributes)
    colnames(node_attributes)[1] = "SYMBOL"
    unique_symbols = unique(node_attributes$SYMBOL)
    node_attributes = node_attributes[node_attributes$SYMBOL %in% unique_symbols, ]  # only select genes that are not duplicated
    node_attributes = rm_row_with_duplicated_value_in_col(node_attributes, "SYMBOL")
    node_attributes = fast_merge(
      nodes,
      node_attributes,
      by.x = "SYMBOL",
      by.y = "SYMBOL",
      all.x = TRUE,
      all.y = FALSE,
      allow.cartesian = FALSE
    )  # this one ensures that the node attributes are assigned to the right nodes.
    #   attr[is.na(attr)] = 0
    #   row.names(attr)= attr$SYMBOL
    attributes_names = colnames(node_attributes)[1:ncol(node_attributes)]
    for (attr in attributes_names)
    {
      graph = set_vertex_attr(graph, name = attr, value = node_attributes[, attr])
    }
  }
  
  
  
  
  #   attr= fast_merge(nodes, attr, by.x="SYMBOL", by.y="SYMBOL")
  #   attr[is.na(attr)] = 0
  #   row.names(attr)= attr$SYMBOL
  
  annotation_sources = get_annotation_sources(excluded_terms = excluded_GO_terms)
  attr_list = list()
  for (i in c(1:length(annotation_sources)))
  {
    attr_list[[i]] = annotation_sources[[i]]$gene_geneset_table
  }
  
  result =  get_network_layout_with_enrichment_score(
    graph,
    clustering_method = clustering_method,
    max_cluster_size = max_cluster_size ,
    ndim = ndim,
    layout_method = layout_method,
    space_ratio = space_ratio,
    node_score_field = node_score_field,
    enrichment_score_field =
      enrichment_score_field,
    sampling_size = sampling_size,
    beta = beta,
    attr = attr_list
  )
  
  
  final_layout = result[[1]]
  final_subnets = result[[2]]
  
  
  cols = c()
  for (col in  graph_attr_names(final_subnets[[1]]))
  {
    if (col != "element_node_names")
    {
      cols = c(cols, col)
    }
  }
  
  
  
  for (subnet in final_subnets)
  {
    subnet_attr = graph_attr(subnet)
    subnet_attr = subnet_attr[cols]
    
    subnet_attr = data.frame(subnet_attr)
    
    
    
    
    if (exists("subnets_attr"))
    {
      subnets_attr = rbind(subnets_attr, subnet_attr)
    } else
    {
      subnets_attr = subnet_attr
    }
  }
  
  if (!is.null(node_score_field))
  {
    subnet_enrichment_df = result[[3]]
    subnets_attr = fast_merge(subnets_attr,
                              subnet_enrichment_df,
                              by.x = "graph_name",
                              by.y = "subnet_name")
    annotation_score_list = result[[4]]
    combined_annotations = combined_word_annotations(
      annotation_score_list,
      annotation_sources,
      subnets = final_subnets,
      n_top_annotation = n_top_annotation,
      importance_field = importance_field,
      n_top_genes = n_top_genes,
      graph = graph
    )
    cl_top_annotation = combined_annotations[[2]]
    cl_top_annotation = cl_top_annotation[, c(1, 2, 3, 5, 6)]  # remove the level column before merging as the subnets_attr already contains this column
    subnets_attr = fast_merge(subnets_attr,
                              cl_top_annotation,
                              by.x = "graph_name",
                              by.y = "subnet_name")
  }
  subnets_attr$enrichment_subnet_score_with_sign  = subnets_attr$enrichment_subnet_score * subnets_attr$change_direction
  
  updated_subnets = list()
  
  for (subnet in final_subnets)
  {
    #print(graph_attr(subnet))
    net_type = graph_attr(subnet)[["net_type"]]
    if (net_type < 2)
      # root and intermediate subnet
    {
      nodes_names = data.frame(name = vertex_attr(subnet)[["name"]])
      
      nodes_attr = fast_merge(nodes_names,
                              subnets_attr,
                              by.x = "name",
                              by.y = "graph_name")
      subnet_name = graph_attr(subnet)[["graph_name"]]
      for (col in colnames(nodes_attr))
      {
        if (col != "n_genes.x" & col != "n_genes.y")
        {
          if (col == "name")
          {
            subnet = set_vertex_attr(subnet, name = "subnet_name", value = subnet_name)
          }
          else
          {
            subnet = set_vertex_attr(subnet, name = col, value =  nodes_attr[, col])
          }
        } else
          if (col == "n_genes.x")
          {
            subnet = set_vertex_attr(subnet, name = "n_genes", value =  nodes_attr[, col])
          }
      }
      
    }
    updated_subnets[[length(updated_subnets) + 1]] = subnet
  }
  
  
  graph_vertices = data.frame(node_name = V(graph)$name)
  graph_vertices_attr = fast_merge(graph_vertices,
                                   final_layout,
                                   by.x = "node_name",
                                   by.y = "node_name")
  graph_vertices_attr$subnet_name = as.character(graph_vertices_attr$subnet_name)
  for (attr in colnames(graph_vertices_attr))
  {
    if (attr != "node_name")
    {
      graph = set_vertex_attr(graph, name = attr, value = graph_vertices_attr[, attr])
    }
  }
  graph = set_graph_attr(graph, "net_type", 0)  # root graph
  # Now concatenating all graphs
  library(plyr)
  attrs = as_data_frame(graph, "vertices")
  el <- as_data_frame(graph)
  for (subnet in updated_subnets)
  {
    #print(graph_attr(subnet)$graph_name)
    if (graph_attr(subnet)$net_type < 2)
    {
      #print("concatenating graph")
      
      attrs_subnet = as_data_frame(subnet, "vertices")
      #attrs_subnet = data.frame(vertex_attr(subnet))
      
      el_subnet = as_data_frame(subnet)
      #print(648)
      attrs = rbind.fill(attrs, attrs_subnet)
      el = rbind.fill(el, el_subnet)
    }
    # else
    # {
    #   #print("final subnet, pass")
    # }
  }
  el = unique(el)
  attrs = unique(attrs)
  attrs_name = data.frame(name = attrs$name)
  attrs = fast_merge(attrs_name, attrs, by.x = "name", by.y = "name")
  graph = graph_from_data_frame(el, directed = FALSE, vertices = attrs)
  if (ndim == 2)
  {
    graph = add_vertices(
      graph,
      nv = 1,
      name = "1",
      radius = 1,
      enrichment_subnet_score = 0,
      X1 = 0,
      X2 = 0
    )
  } else
  {
    graph = add_vertices(
      graph,
      nv = 1,
      name = "1",
      radius = 1,
      enrichment_subnet_score = 0,
      X1 = 0,
      X2 = 0,
      X3 = 0
    )
  }
  
  
  for (i in c(1:nrow(subnets_attr)))
  {
    split_top_annotation = unlist(strsplit(subnets_attr[i, "top_annotation"], split =
                                             "\n"))
    subnets_attr[i, "top_annotation"] = paste(split_top_annotation, collapse = "; ")
    split_annotation = unlist(strsplit(subnets_attr[i, "annotation"], split =
                                         "\n"))
    subnets_attr[i, "annotation"] = paste(split_annotation, collapse = "; ")
  }
  
  
  sign_subnets_attr = subnets_attr[subnets_attr$subnet_qvalue < 0.05 & !is.na(subnets_attr$subnet_qvalue),]
  if (is.null(final_layout$is_external[1]))
  {
    final_layout = get_external_compartment_field(final_layout)
  }
  
  print("1521 final layout")
  print(final_layout)
  out = list(
    graph = graph,
    layout = final_layout,
    subnets = updated_subnets,
    subnets_attr = subnets_attr,
    annotations = combined_annotations,
    sign_subnets = sign_subnets_attr,
    annotation_sources = annotation_sources
  )
  return(out)
}


# pipeline_out: result from get_network_layout_pipeline
save_layout = function(outdir = "network",
                       gene_description_file = paste0(data_path, "/all_gene_summary.csv"),
                       pipeline_out,
                       selected_subnets = NULL,
                       #c("1.4.3.4", "1.5.5.3"), # immune subnet,
                       selected_genes = NULL,
                       #c("NKX6-2","ASCL1","MYCN","BASP1"),
                       qvalue_cutoff = 0.05,
                       cut_off_field = "logFC",
                       sort_fields = c(cut_off_field, "score"),
                       cut_off_value = 1,
                       selected_fields = c(
                         "node_name",
                         "score",
                         cut_off_field,
                         "description",
                         "summary",
                         "subnet_name",
                         "degree",
                         "pvalue",
                         "fdr",
                         "annotation"
                       ),
                       draw_heatmap = TRUE,
                       gene_expression_table = NULL,
                       #"~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_log.txt",
                       sample_group_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_sample_description.csv",
                       displayed_groups = c("GSC", "NSC", "NHA"),
                       n_top_genes = 50,
                       importance_field = "score",
                       is_log_gene_expression_table = FALSE,
                       do_extract_subnet_data = FALSE,
                       baseline_samples = NULL,
                       use_cluster = FALSE,
                       n_root_master_regulators = 5,
                       save_network_layout_only = FALSE,
                       export_top_master_genes_data = FALSE,
                       debug = TRUE,
                       baseline_table = NULL,
                       layout_dir = "network/layout/1")
{
  #' @param n_root_master_regulators interger number of top level master regulators to extract data for each master regulator
  # unlist to get input data
  graph = pipeline_out[[1]]
  final_layout = pipeline_out[[2]]
  subnets = pipeline_out[[3]]
  subnets_attr = pipeline_out[[4]]
  annotations = pipeline_out[[5]]
  sign_subnets_attr = pipeline_out[[6]]
  annotation_sources = pipeline_out[[7]]
  
  combined_annotations = annotations[[1]]
  subnet_top_annotation = annotations[[2]]
  subnet_annotation_by_coverage = annotations[[3]]
  combined_subnet_qvalues = annotations[[4]]
  combined_gene_geneset_table = get_combined_gene_geneset_table(annotation_sources)
  
  dir.create(outdir)
  dir.create(layout_dir)
  if (debug == TRUE)
  {
    
    inputs = list(outdir = outdir,
     gene_description_file = gene_description_file,
     pipeline_out = pipeline_out,
     selected_subnets =  selected_subnets,
     #c("1.4.3.4", "1.5.5.3"), # immune subnet,
     selected_genes = selected_genes,
     #c("NKX6-2","ASCL1","MYCN","BASP1"),
     qvalue_cutoff =qvalue_cutoff,
     cut_off_field = cut_off_field,
     sort_fields = sort_fields,
     cut_off_value = cut_off_value,
     selected_fields = selected_fields,
     draw_heatmap = draw_heatmap,
     gene_expression_table = gene_expression_table, 
     #"~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_log.txt",
     sample_group_table = sample_group_table ,
     displayed_groups = displayed_groups,
     n_top_genes = n_top_genes,
     importance_field = importance_field,
     is_log_gene_expression_table = is_log_gene_expression_table,
     do_extract_subnet_data = do_extract_subnet_data,
     baseline_samples = baseline_samples,
     use_cluster = use_cluster,
     n_root_master_regulators = n_root_master_regulators,
     save_network_layout_only = save_network_layout_only,
     export_top_master_genes_data =  export_top_master_genes_data ,
     debug = debug,
     baseline_table = baseline_table,
     layout_dir=layout_dir)
    
    
    outfile = paste(layout_dir, "/save_layout.RDS", sep = "")
    saveRDS(inputs, file = outfile)
  }
  graph_file = paste(layout_dir, "/out_graph.graphml", sep = "")
  write_graph(graph, file = graph_file, format = "graphml")
  if (!save_network_layout_only)
  {
    get_hierarchy_network(
      final_layout = final_layout,
      subnets_attr = subnets_attr,
      outdir = outdir,
      diff_exp_field = cut_off_field
    )
    
    final_layout = get_gene_description(
      final_layout,
      subnets_attr,
      gene_description_file = paste0(data_path, "/all_gene_summary.csv")
    )
    
    
    
    annotation_dir = paste(outdir, "/annotations_data", sep = "")
    dir.create(annotation_dir)
    
    
    for (i in c(1:length(combined_annotations)))
    {
      annontation_subnet_table = combined_annotations[[i]]
      sn_an_final = annontation_subnet_table[[1]]
      sn_an_sum =  annontation_subnet_table[[2]]
      sn_an_by_coverage =  annontation_subnet_table[[3]]
      
      filename = paste(annotation_dir, "/subnet_annotation_", i, ".txt", sep =
                         "")
      write.table(
        sn_an_final,
        filename,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      #filename = paste(outdir, "/subnet_annotation_sum_", i,".txt", sep="")
      #write.table( sn_an_sum, filename, sep="\t", row.names=FALSE, quote=FALSE)
      filename = paste(annotation_dir,
                       "/subnet_annotation_by_coverage_",
                       i,
                       ".txt",
                       sep = "")
      write.table(
        sn_an_by_coverage ,
        filename,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      
    }
    
    filename = paste(annotation_dir, "/subnet_top_annotation_.txt", sep =
                       "")
    write.table(
      subnet_top_annotation,
      filename,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    filename = paste(annotation_dir,
                     "/combined_subnet_annotation_by_coverage.txt",
                     sep = "")
    write.table(
      subnet_annotation_by_coverage,
      filename,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    filename = paste(annotation_dir,
                     "/combined_subnet_qvalues_annotation.txt",
                     sep = "")
    write.table(
      combined_subnet_qvalues,
      filename,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    
    
    
    subnets_attr_file = paste(outdir, "/subnets_attr.txt", sep = "")
    write.table(
      subnets_attr,
      file = subnets_attr_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    layout_file = paste(outdir, "/layout.txt", sep = "")
    write.table(
      final_layout,
      file = layout_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    sign_subnets_file = paste(outdir, "/sign_subnets.txt", sep = "")
    write.table(
      sign_subnets_attr,
      file = sign_subnets_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    
    external_compartment = external_compartment_analysis(final_layout,
                                                         cut_off_field = cut_off_field,
                                                         cut_off_value = cut_off_value)
    external_compartment_file = paste(outdir, "/external_compartment.txt", sep =
                                        "")
    write.table(
      external_compartment[[1]],
      file = external_compartment_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    external_compartment_short_file = paste(outdir, "/external_compartment_short.txt", sep =
                                              "")
    external_compartment_short = external_compartment[[2]]
    write.table(
      external_compartment_short,
      file = external_compartment_short_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    
    
    sign_subnets_attr_l3 =  sign_subnets_attr[sign_subnets_attr$level <= 3, c(1, 2, 3, 5, 12, 15, 17)]
    colnames(sign_subnets_attr_l3) = c(
      "Subnet",
      "n_genes",
      "level",
      "Increase/Decrease",
      "qvalue",
      "annotation",
      "Masters"
    )
    
    sign_subnets_attr_l3_file = paste(outdir, "/sign_subnets_l3.txt", sep =
                                        "")
    write.table(
      sign_subnets_attr_l3,
      file = sign_subnets_attr_l3_file ,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    if (!is.null(selected_subnets))
    {
      for (selected_subnet in selected_subnets)
      {
        tryCatch({
          extract_subnet_data(
            final_layout = final_layout,
            subnet = selected_subnet,
            subnets_attr = subnets_attr,
            gene_description_file = gene_description_file,
            qvalue_cutoff = 0.05,
            outdir = paste(outdir, "/",selected_subnet, sep =
                             ""),
            sort_fields = sort_fields,
            cut_off_field = cut_off_field,
            cut_off_value = cut_off_value,
            selected_fields = selected_fields,
            subnet_annotation_by_coverage = subnet_annotation_by_coverage,
            combined_subnet_qvalues = combined_subnet_qvalues,
            combined_gene_geneset_table = combined_gene_geneset_table,
            draw_heatmap = draw_heatmap,
            gene_expression_table = gene_expression_table,
            sample_group_table = sample_group_table,
            displayed_groups = displayed_groups,
            n_top_genes = n_top_genes,
            importance_field = importance_field,
            is_log_gene_expression_table = is_log_gene_expression_table,
            baseline_samples = baseline_samples,
            use_cluster = use_cluster,
            annotation_sources = annotation_sources,
            baseline_table = baseline_table
          )
        })
      }
    }
    
    sorted_layout = final_layout[order(final_layout$score, decreasing = TRUE),]
    if (export_top_master_genes_data == TRUE)
    {
      top_master_genes = sorted_layout$node_name[1:n_root_master_regulators]
      
      for (gene in top_master_genes)
      {
        extract_gene_function(
          gene_name = gene,
          layout = final_layout,
          subnets_attr = subnets_attr,
          outdir = paste(outdir, "/Masters", sep =
                           ""),
          subnet_annotation_by_coverage = subnet_annotation_by_coverage,
          combined_subnet_qvalues = combined_subnet_qvalues,
          combined_gene_geneset_table = combined_gene_geneset_table,
          diff_exp_field = cut_off_field,
          selected_fields = selected_fields,
          draw_heatmap = draw_heatmap,
          gene_expression_table = gene_expression_table,
          sample_group_table = sample_group_table,
          displayed_groups = displayed_groups,
          n_top_genes = n_top_genes,
          importance_field = importance_field,
          is_log_gene_expression_table = is_log_gene_expression_table,
          n_top_pathway_to_draw_heatmap = 5,
          baseline_samples = baseline_samples,
          annotation_sources = annotation_sources,
          use_cluster = use_cluster,
          baseline_table = baseline_table
        )
      }
      
    }
    
    
    if (!is.null(selected_genes))
    {
      for (gene in selected_genes)
      {
        extract_gene_function(
          gene_name = gene,
          layout = final_layout,
          subnets_attr = subnets_attr,
          outdir = paste(outdir, "/selected_genes", sep =
                           ""),
          subnet_annotation_by_coverage = subnet_annotation_by_coverage,
          combined_subnet_qvalues = combined_subnet_qvalues,
          combined_gene_geneset_table = combined_gene_geneset_table,
          diff_exp_field = cut_off_field,
          selected_fields = selected_fields,
          draw_heatmap = draw_heatmap,
          gene_expression_table = gene_expression_table,
          sample_group_table = sample_group_table,
          displayed_groups = displayed_groups,
          n_top_genes = n_top_genes,
          importance_field = importance_field,
          is_log_gene_expression_table = is_log_gene_expression_table,
          n_top_pathway_to_draw_heatmap = 5,
          baseline_samples = baseline_samples,
          annotation_sources = annotation_sources,
          use_cluster = use_cluster,
          baseline_table = baseline_table
        )
      }
    }
    
    
    
    
    
    
    
    
    subnets_l123 = as.character(sign_subnets_attr_l3$Subnet)
    if (do_extract_subnet_data == TRUE)
    {
      subnets_l123 = c("1", subnets_l123)
    }
    else {
      subnets_l123 = c("1")
    }
    subnets_dir = paste(outdir, "/subnets", sep = "")
    dir.create(subnets_dir)
    for (sub in subnets_l123)
    {
      print(paste("extracting subnet data:", sub))
      extract_subnet_data(
        final_layout = final_layout,
        subnet = sub,
        subnets_attr = subnets_attr,
        gene_description_file = gene_description_file,
        qvalue_cutoff = 0.05,
        outdir = paste(subnets_dir,"/",sub, sep = ""),
        sort_fields = sort_fields,
        cut_off_field = cut_off_field,
        cut_off_value = cut_off_value,
        selected_fields = selected_fields,
        subnet_annotation_by_coverage = subnet_annotation_by_coverage,
        combined_subnet_qvalues = combined_subnet_qvalues,
        combined_gene_geneset_table = combined_gene_geneset_table,
        draw_heatmap = draw_heatmap,
        gene_expression_table = gene_expression_table,
        sample_group_table = sample_group_table,
        displayed_groups = displayed_groups,
        n_top_genes = n_top_genes,
        importance_field = importance_field,
        baseline_samples = baseline_samples,
        use_cluster = use_cluster,
        annotation_sources = annotation_sources,
        baseline_table = baseline_table
      )
    }
    
    #print("finished loop")
    
    
    
    out = list(
      graph = graph,
      layout = final_layout,
      subnets = subnets,
      subnets_attr = subnets_attr,
      annotations = annotations,
      sign_subnets = sign_subnets_attr,
      annotation_sources = annotation_sources,
      combined_annotations = combined_annotations,
      subnet_top_annotation = subnet_top_annotation,
      subnet_annotation_by_coverage = subnet_annotation_by_coverage,
      combined_subnet_qvalues = combined_subnet_qvalues,
      combined_gene_geneset_table = combined_gene_geneset_table
    )
    
    outfile = paste(layout_dir,
                    "/get_network_layout_pipeline_after_saving.RDS",
                    sep = "")
    saveRDS(out, file = outfile)
    return(out)
  }
}


# from annotation sources, extract gene_geneset tables and concatenate columwise
get_combined_gene_geneset_table = function(annotation_sources)
{
  for (i in c(1:length(annotation_sources)))
  {
    gene_geneset_table = annotation_sources[[i]]$gene_geneset_table
    if (!exists("combined_gene_geneset_table"))
    {
      combined_gene_geneset_table = gene_geneset_table
    } else
    {
      combined_gene_geneset_table = fast_merge(
        combined_gene_geneset_table,
        gene_geneset_table,
        by.x = "SYMBOL",
        by.y = "SYMBOL",
        all.x = TRUE,
        all.y = TRUE
      )
    }
  }
  return(combined_gene_geneset_table)
  
}

# add gene name, description, subnet attributes fields to final layout
get_gene_description = function(final_layout,
                                subnets_attr,
                                gene_description_file = paste0(data_path, "/all_gene_summary.csv"))
{
  if (!file.exists(gene_description_file))
  {
    get_all_gene_info(organism = "Homo sapiens",
                      batch_size = 300,
                      outfile = gene_description_file)
  }
  
  library(readr)
  all_gene_summary <-
    read_delim(gene_description_file,
               "\t",
               escape_double = FALSE,
               trim_ws = TRUE)
  all_gene_summary = unique(all_gene_summary[, c(2, 5, 6)])
  final_layout = fast_merge(
    final_layout,
    all_gene_summary,
    by.x = "node_name",
    by.y = "SYMBOL",
    all.x = TRUE
  )
  final_layout = fast_merge(
    final_layout,
    subnets_attr[, c("graph_name", "annotation")],
    by.x = "subnet_name",
    by.y = "graph_name",
    all.x = TRUE
  )
  return(final_layout)
}



get_network_layout_with_enrichment_score = function(graph,
                                                    clustering_method = "louvain",
                                                    max_cluster_size = 100 ,
                                                    ndim = 3,
                                                    layout_method = "fr",
                                                    space_ratio,
                                                    node_score_field = "logFC",
                                                    enrichment_score_field =
                                                      NULL,
                                                    sampling_size = 1000,
                                                    beta = 0.5,
                                                    attr = NULL)
{
  result = get_network_layout(
    graph = graph,
    clustering_method = clustering_method,
    max_cluster_size = max_cluster_size,
    ndim = ndim,
    layout_method = layout_method,
    space_ratio = space_ratio,
    node_score_field = node_score_field,
    beta = beta,
    attr = attr
  )
  
  final_layout = result[[1]]
  final_subnets = result[[2]]
  
  node_attr = data.frame(vertex_attr(graph))
  
  node_attr[is.na(node_attr)] = 0
  
  final_layout = fast_merge(final_layout, node_attr, by.x = "node_name", by.y =
                              "name")
  if (!is.null(node_score_field))
  {
    if (!is.null(enrichment_score_field))
    {
      for (field in enrichment_score_field)
      {
        final_layout[, paste("enrichment", field, sep = "_")] = get_enrichment_score(final_layout[, field], sampling_size = sampling_size)
      }
    }
    node_score_field = "subnet_score"
    subnet_enrichment_df = get_subnet_enrichment_score_from_subnets_list(final_subnets,
                                                                         node_attrs = final_layout[, c("node_name", node_score_field)],
                                                                         sampling_size = sampling_size)
    subnet_pvalues =  subnet_enrichment_df[[1]]
    subnet_qvalues =   subnet_enrichment_df[[2]]
    subnet_enrichment_values =  subnet_enrichment_df[[3]]
    subnet_enrichment_df = data.frame(
      subnet_name = subnet_pvalues$subnet_name,
      subnet_pvalue = subnet_pvalues$subnet_score,
      subnet_qvalue = subnet_qvalues$subnet_score,
      enrichment_subnet_score = subnet_enrichment_values$subnet_score
    )
    #final_layout = get_subnet_enrichment_score(final_layout, sampling_size = sampling_size)
    final_layout = fast_merge(final_layout,
                              subnet_enrichment_df,
                              by.x = "subnet_name",
                              by.y = "subnet_name")
    final_layout$enrichment_subnet_score_with_sign =  abs(final_layout$enrichment_subnet_score) * sign(final_layout[, node_score_field])  # something wrong here
    
  }
  
  if (!is.null(attr))
  {
    node_attr_score_df = result[[3]]
    return(list(
      final_layout,
      final_subnets,
      subnet_enrichment_df,
      node_attr_score_df
    ))
  } else
  {
    return(list(final_layout, final_subnets, subnet_enrichment_df))
  }
}



# calculate enrichment score for each subnet in subnets,
# subnets: list of all subnets in the graph, node_attrs: data frame of node attributes, first column is the gene name, other columns are attributes and their values
#return 3 data frame , subnet_pvalue, subnet_qvalue, subnet_enrichment_value
# the first column is subnet name, second is the number of genes in subnet,  other columns are attributes (GO pathway)
get_subnet_enrichment_score_from_subnets_list = function(subnets,
                                                         node_attrs,
                                                         sampling_size = 1000,
                                                         use_random_sampling = FALSE)
{
  colnames(node_attrs)[1] = "node_name"  # change the name of first column of node_attrs
  n_genes_list = c()
  
  subnet_names = c()
  subnet_levels = c()
  n_subnet = length(subnets)
  for (subnet in subnets)
  {
    subnet_name =  graph_attr(subnet)[["graph_name"]]
    subnet_names = c(subnet_names, subnet_name)
    n_genes = graph_attr(subnet)[["n_genes"]]
    n_genes_list = c(n_genes_list, n_genes)
    subnet_level =  graph_attr(subnet)[["level"]]
    subnet_levels = c(subnet_levels, subnet_level)
  }
  subnet_pvalues = data.frame(subnet_name = subnet_names,
                              level = subnet_levels,
                              n_genes = n_genes_list)
  subnet_qvalues = data.frame(subnet_name = subnet_names,
                              level = subnet_levels,
                              n_genes = n_genes_list)
  subnet_enrichment_values = data.frame(subnet_name = subnet_names,
                                        level = subnet_levels,
                                        n_genes = n_genes_list)
  rownames(subnet_pvalues) = subnet_names
  rownames(subnet_qvalues) = subnet_names
  rownames(subnet_enrichment_values) = subnet_names
  i = 0
  time_zero = Sys.time()
  n_col = ncol(node_attrs)
  cols = colnames(node_attrs)
  for (col_i in c(2:n_col))
    # process each attribute, except the fisrt column that is gene name
  {
    col = cols[col_i]
    i = i + 1
    start_time = Sys.time()
    attr_scores = abs(node_attrs[, col_i])
    attr_scores[is.na(attr_scores)] = 0.
    attr_pvalues = vector(length = n_subnet)
    #subnet_names = c()
    
    M = mean(attr_scores)
    S = sd(attr_scores)
    node_attrs_attr = node_attrs[, col_i]
    for (subnet_i in c(1:n_subnet))
    {
      subnet = subnets[[subnet_i]]
      element_nodes = graph_attr(subnet, name = "element_node_names")
      match_index = fmatch(element_nodes, node_attrs$node_name)
      subnet_node_scores = node_attrs_attr[match_index]
      
      if (use_random_sampling)
      {
        random_data = c()
        
        for (i in c(1:sampling_size))
        {
          random_subnet_score = mean(sample(attr_scores, size = n, replace = FALSE))
          random_data = c(random_data, random_subnet_score)
        }
        
        p_value = pnorm(
          subnet_score,
          mean = mean(random_data),
          sd = sd(random_data),
          lower.tail = FALSE
        )
      } else
      {
        m = mean(subnet_node_scores)  # subnet mean
        s = sd(subnet_node_scores)  # subnet sd
        n = graph_attr(subnet)[["n_genes"]]
        t_value =  (m - M) / sqrt(s ^ 2 / n + S ^ 2 / n)   # GAGE two -sample t -test
        df = (n - 1) * (s ^ 2 + S ^ 2) ^ 2 / (s ^ 4 + S ^ 4)
        p_value = pt(t_value, df = df, lower.tail = FALSE)
      }
      
      attr_pvalues[subnet_i] = p_value
    }
    attr_qvalues = p.adjust(attr_pvalues, method = "BH")
    
    attr_qvalues = attr_qvalues + 1e-6
    attr_qvalues[attr_qvalues > 1] = 1
    attr_enrichment_values = -log10(attr_qvalues) / 6
    subnet_pvalues[, col] = attr_pvalues
    subnet_qvalues[, col] = attr_qvalues
    subnet_enrichment_values[, col] = attr_enrichment_values
    end_time = Sys.time()
    elapsed = end_time - start_time
    total_elapsed = end_time - time_zero
    # print(
    #   paste(
    #     i,
    #     " done calculating attribute enrichment score of",
    #     col,
    #     "for",
    #     elapsed,
    #     "Total elapsed time:",
    #     total_elapsed
    #   )
    # )
  }
  
  
  out = list(subnet_pvalues, subnet_qvalues, subnet_enrichment_values)
  #out = data.frame(subnet_name=subnet_names, n_genes =n_genes_list, subnet_pvalue=subnet_pvalues, subnet_qvalue =subnet_qvalues, subnet_enrichment_value = subnet_enrichment_values)
  return(out)
}

# Input is a dataframe with column (subnet_name,  n_genes, attr1, attr2, etc.) where attr1 is the qvalue of this attr for each subnet
#For each attributes, select subnets that is significant. Among those significant subnets, assign attributes to a subnet or several subnets that meet criteria:
#
subnet_annotation = function(subnet_qvalues,
                             sign_value = 0.05,
                             max_split = 3)
{
  attrs = colnames(subnet_qvalues)[4:ncol(subnet_qvalues)]
  rownames(subnet_qvalues) = subnet_qvalues$subnet_name
  subnet_ann = subnet_qvalues
  subnet_ann[, c(4:ncol(subnet_ann))] = 0
  for (attr in attrs)
  {
    subnet_attr_ann = subnet_qvalues[subnet_qvalues[, attr] < sign_value, c(1:3)]
    subnet_attr_ann_attr = subnet_qvalues[subnet_qvalues[, attr] < sign_value, attr]
    subnet_attr_ann = cbind(subnet_attr_ann, attr = subnet_qvalues[subnet_qvalues[, attr] < sign_value, attr])
    
    levels_count = ddply(subnet_attr_ann, ~ level,  nrow)
    levels_count = levels_count[order(levels_count$level), ]
    row_i = 1
    level = NULL
    # The level to assign the attribute is the highest level one still has the number of  significant pathways equal 1. And this pathway at this level is assigned to the annotation
    while (row_i <= nrow(levels_count))
    {
      count = levels_count[row_i, "V1"]
      if (count == 1)
      {
        if (row_i < nrow(levels_count))
        {
          row_i = row_i + 1
        } else
        {
          level = levels_count[row_i, "level"] # reach the bottom level
          break()
        }
      }
      else
      {
        level = levels_count[row_i, "level"]
        if (levels_count[row_i, "V1"] > max_split)
          # Maximum number of subnets at this level allowed
        {
          level = level - 1
        }
        if (level == 1)
        {
          level = NULL
        }
        break()
      }
      
    }
    
    if (!is.null(level))
    {
      combined = TRUE
      if (combined)
      {
        if (level > 2)
        {
          levels = c(level, level - 1)
        } else
        {
          levels = c(level)
        }
      } else
      {
        levels = c(level)
      }
      
      subnet_attr_ann = subnet_attr_ann[subnet_attr_ann$level %in%  levels, ]
      #print(paste(attr, "accepted levels:"))
      #print(levels)
      #print(levels_count)
      subnet_ann[subnet_ann$subnet_name %in% subnet_attr_ann$subnet_name, attr] = 1
    }
  }
  return(subnet_ann)
}

# Input is a dataframe with column (subnet_name,  n_genes, attr1, attr2, etc.) where attr1 is the qvalue of this attr for each subnet
#For each attributes, if number of sign subnet at a level =1, continue until split. From split point, select the subnet from the level above, and for each level select one subnet that has lowest qvalue and is significant.
#So every level will have one subnet selected from split point for each attribute
subnet_annotation_level_by_level = function(subnet_qvalues,
                                            sign_value = 0.05,
                                            subnets,
                                            gene_geneset_table)
{
  attrs = colnames(subnet_qvalues)[4:ncol(subnet_qvalues)]  # list of attribute names
  subnet_qvalues$subnet_name = as.character(subnet_qvalues$subnet_name)
  rownames(subnet_qvalues) = subnet_qvalues$subnet_name
  subnet_ann = subnet_qvalues
  subnet_ann = subnet_ann[subnet_ann$level > 1, ]
  subnet_ann[, c(4:ncol(subnet_ann))] = 0
  node_subnet_dict = new.env(hash = TRUE, size = length(subnets))
  for (subnet in subnets)
  {
    subnet_name = graph_attr(subnet, name = "graph_name")
    subnet_genes = graph_attr(subnet, name = "element_node_names")
    node_subnet_dict[[subnet_name]] = subnet_genes
  }
  for (attr in attrs)
  {
    #attr = attrs[1]
    subnet_attr_ann = subnet_qvalues[subnet_qvalues[, attr] < sign_value, c(1:3)]
    subnet_attr_ann_attr = subnet_qvalues[subnet_qvalues[, attr] < sign_value, attr]
    subnet_attr_ann = cbind(subnet_attr_ann, attr = subnet_qvalues[subnet_qvalues[, attr] < sign_value, attr])
    
    levels_count = ddply(subnet_attr_ann, ~ level,  nrow)
    levels_count = levels_count[order(levels_count$level), ]
    row_i = 1
    level = NULL
    # Determine the level one step above split point
    while (row_i <= nrow(levels_count))
    {
      count = levels_count[row_i, "V1"]
      if (count == 1)
      {
        if (row_i < nrow(levels_count))
        {
          row_i = row_i + 1
        } else
        {
          level = levels_count[row_i, "level"] # reach the bottom level
          break()
        }
      }
      else
      {
        level = levels_count[row_i, "level"]  # reaching slit point
        #if ( levels_count[row_i,"V1"] > max_split)  # Maximum number of subnets at this level allowed
        #{
        level = level - 1
        #levels = c(level, level +1)
        #}
        #         if (level == 1)
        #         {
        #           level = NULL # if split begin at root, then the attribute is too broad, skip
        #           # level = 2 # if split begin at root, then use lower level
        #           #level = 1
        #           #levels = c(1,2)
        #         }
        break()
      }
      
    }
    
    levels = levels_count[levels_count$level >= level, "level"]
    #levels = c(level, level + 1)
    #levels = c(level)
    if (is.null(level))
    {
      print(levels_count)
      next
    }
    #print(attr)
    parent_clusters = NULL
    attr_genes = gene_geneset_table[gene_geneset_table[, attr] == 1, 1]
    for (lev in levels)
    {
      subnet_attr_ann_level = subnet_attr_ann[subnet_attr_ann$level == lev, ]
      if (!is.null(parent_clusters))
        # select only the subnets that belong to the parent subnets
      {
        subnet_attr_ann_level = subnet_attr_ann_level[substr(
          subnet_attr_ann_level$subnet_name,
          1,
          nchar(subnet_attr_ann_level$subnet_name) - 2
        ) %in% parent_clusters,]
      }
      
      subnet_names = subnet_attr_ann_level$subnet_name
      max_coverage = 0
      selected_subnet_name = NULL
      for (subnet_name in subnet_names)
      {
        subnet_genes = node_subnet_dict[[subnet_name]]
        n_subnet_genes = length(subnet_genes)
        n_overlapped_genes = number_of_overlapped_genes(attr_genes, subnet_genes)
        coverage = n_overlapped_genes / n_subnet_genes
        if (coverage >=  max_coverage)
        {
          max_coverage = coverage
          selected_subnet_name = subnet_name
        }
      }
      
      subnet_ann[subnet_ann$subnet_name == selected_subnet_name, attr] = 1
      #print(paste(attr, lev))
      #print(selected_subnet_name)
      parent_clusters = c(selected_subnet_name)
    }
    
  }
  sn_an_sum = data.frame(n_clusters = colSums(subnet_ann[4:ncol(subnet_ann)]))
  sn_an_sum$attribute_name = rownames((sn_an_sum))
  sn_an_sum = sn_an_sum[sn_an_sum$n_clusters > 0, ]
  selected_attributes = sn_an_sum$attribute_name
  subnet_ann = subnet_ann[, selected_attributes]
  subnet_ann = cbind(data.frame(subnet_name = rownames(subnet_ann)), subnet_ann)
  return(subnet_ann)
}

# fast count the nunmber of overlapped genes in two geneset
number_of_overlapped_genes = function(geneset1, geneset2)
{
  n_share_genes = sum(geneset1 %fin% geneset2)
  return(n_share_genes)
}

# select annotation attributes that contain more than min_ngene
select_attr_by_number_gene = function(attr, min_ngene = 5)
{
  attr_ngenes = data.frame(n_gene = colSums(attr[2:ncol(attr)]))
  attr_ngenes$attr_id = rownames((attr_ngenes))
  # only select attributes that contain more than min_ngene
  attr_ngenes = attr_ngenes[attr_ngenes$n_gene >= min_ngene, ]
  attr_cols = attr_ngenes$attr_id
  attr = cbind(SYMBOL = attr[, 1], attr[, attr_cols])
  attr$SYMBOL = as.character(attr$SYMBOL)
  return(attr)
}

# subnet is assigned to an attribute if it's the lowest level subnet that contain more than predefined coverage fraction(0.5) of number of genes in this attribute
get_subnet_annotation_by_coverage = function(subnets,
                                             attr,
                                             coverage = 0.5,
                                             min_ngene = 5)
{
  n_attr = ncol(attr) - 1
  subnet_names = c()
  for (subnet in subnets)
  {
    subnet_names = c(subnet_names, graph_attr(subnet, name = "graph_name"))
    
  }
  sub_an_df = data.frame(subnet = subnet_names)
  n_subnet = nrow(sub_an_df)
  sub_an_df[, attr_cols] = 0
  rownames(sub_an_df) = subnet_names
  for (attr_i in c(1:n_attr))
  {
    at = attr[, c(1, attr_i + 1)]
    at_ngene = attr_ngenes[attr_i, "n_gene"]
    for (subnet_i in c(1:n_subnet))
    {
      subnet = subnets[[subnet_i]]
      subnet_genes = graph_attr(subnet, name = "element_node_names")
      subnet_genes_match = fmatch(subnet_genes, at$SYMBOL)
      subnet_at = at[subnet_genes_match, 2]
      coverage_score = sum(subnet_at) / at_ngene
      
      if (coverage_score > coverage)
      {
        #print(coverage_score)
        sub_an_df[subnet_i, attr_i + 1] = 1
      }
      
    }
    
  }
  # table of the number of subnets that contain more than coverage of particular annotation
  sub_an_df_sum = data.frame(count = colSums(sub_an_df[2:ncol(sub_an_df)]))
  sub_an_df_sum$subnet = rownames(sub_an_df_sum)
  # eliminate annotation that contain only one subnet (root subnet)
  sub_an_df_sum = sub_an_df_sum[sub_an_df_sum$count > 1, ]
  # selectec annoation
  selected_atts = sub_an_df_sum$subnet
  sub_an_df = sub_an_df[, selected_atts]
  sub_an_df_t = t(sub_an_df)
  sub_an_df_t_sum = data.frame(n_go = colSums((sub_an_df_t)))
}


# from result of subnet_annotation, for each attributes, if number of pathways >1 then the attributs is not specific,
#and set up the assignment values for all pathway equal zero for this attributes

trim_subnet_annotation = function(subnet_annotation_df)
{
  attrs = colnames(subnet_qvalues)[4:ncol(subnet_qvalues)]
  rownames(subnet_qvalues) = subnet_qvalues$subnet_name
  for (attr in attrs)
  {
    if (sum(subnet_annotation_df[, attr] > 1))
    {
      subnet_annotation_df[, attr] = 0
    }
  }
  
  return(subnet_annotation_df)
}



plot_network = function(network_layout, cex, scale_ratio = 1)
{
  ndim = ncol(network_layout) - 2
  cex = network_layout$degree * cex
  if (ndim == 2)
  {
    plot(network_layout$X1 * scale_ratio,
         network_layout$X2 * scale_ratio,
         cex = cex)
    text(
      network_layout$X1 * scale_ratio + cex,
      network_layout$X2 * scale_ratio + cex ,
      labels = network_layout$node_name,
      cex = cex
    )
  } else
  {
    if (ndim == 3)
    {
      library("car")
      p = scatter3d(
        x = network_layout$X1 * scale_ratio,
        y = network_layout$X2 * scale_ratio,
        z = network_layout$X3 * scale_ratio,
        point.col = "blue",
        surface = FALSE,
        grid = FALSE,
        labels = network_layout$node_name
      )
      #library(plotly)
      #p = plot_ly(network_layout, x = ~X1, y =~X2, z = ~X3)
      #add_markers(p, symbol = ~node_name)
      #print(p)
    }
  }
}


find_beta_network_layout = function(network_file = "nha_gsc_net.graphml" ,
                                    graph_format = "graphml",
                                    clustering_method = "louvain",
                                    max_cluster_size = 100 ,
                                    ndim = 3,
                                    layout_method = "fr",
                                    space_ratio,
                                    node_score_field = "logFC",
                                    enrichment_score_field = c("node_smooth_score"),
                                    sampling_size = 1000,
                                    n_beta_interval = 11)
{
  print("reading graph file ...")
  graph = read_graph(file = network_file, format = graph_format)
  print("graph file has been read")
  graph = as.undirected(graph, mode = c("collapse"))
  
  for (beta in seq(from = 0,
                   to = 1,
                   length.out = n_beta_interval))
  {
    final_layout =  get_network_layout_with_enrichment_score(
      graph,
      clustering_method = clustering_method,
      max_cluster_size = max_cluster_size ,
      ndim = ndim,
      layout_method = layout_method,
      space_ratio =
        space_ratio,
      node_score_field = node_score_field,
      enrichment_score_field =
        enrichment_score_field,
      sampling_size = sampling_size,
      beta = beta
    )
    enrichment_subnet_score = final_layout[, c("node_name", "enrichment_subnet_score")]
    colnames(enrichment_subnet_score)[2] = as.character(beta)
    if (exists("end_layout"))
    {
      end_layout = fast_merge(end_layout,
                              enrichment_subnet_score,
                              by.x = "node_name",
                              by.y = "node_name")
    } else
    {
      end_layout = enrichment_subnet_score
    }
    
    #print(paste("drawing plot for beta", beta))
    #print(colnames(end_layout))
    if (beta == 0)
    {
      plot(density(end_layout[, as.character(beta)]), xlim = range(-0.2, 1.2))
    } else
    {
      lines(density(end_layout[, as.character(beta)]))
    }
  }
  
  return(end_layout)
}

get_human_go = function(min_ngene = 5,
                        max_ngene = 500,
                        outprefix = paste0(data_path, "/GO"))
{
  library(AnnotationHub)
  ah = AnnotationHub()
  orgs = subset(ah, ah$rdataclass == "OrgDb")
  human = query(orgs, "Homo sapiens")[[1]]
  library(GO.db)
  go = select(
    human,
    keys = keys(human, keytype = "SYMBOL"),
    columns = c("SYMBOL", "ONTOLOGY", "GO"),
    keytype = "SYMBOL"
  )
  go = go[!is.na(go$GO), ]
  goid_key = keys(GO.db, keytype = "GOID")
  go_def = select(
    GO.db,
    keys = goid_key,
    columns = c("DEFINITION", "GOID", "ONTOLOGY", "TERM"),
    keytype = "GOID"
  )
  go = fast_merge(go, go_def[, c("GOID", "TERM")], by.x = "GO", by.y = "GOID")
  human_GO_list = unique(go$GO)
  human_go_def = go_def[go_def$GOID %in% human_GO_list, ]
  
  go_unique = unique(go[, c("GO", "SYMBOL")])
  genes_GO = data.frame(SYMBOL = unique(go_unique$SYMBOL))
  rownames(genes_GO) = genes_GO$SYMBOL
  for (i in c(1:nrow(go_unique)))
  {
    #print(paste("processing row", i))
    SYMBOL = go_unique[i, "SYMBOL"]
    GO = go_unique[i, "GO"]
    genes_GO[SYMBOL, GO] = 1
  }
  genes_GO[is.na(genes_GO)] = 0
  
  filtered_result = filter_annotation_by_ngene(
    genes_GO,
    human_go_def,
    min_ngene = min_ngene,
    max_ngene = max_ngene,
    outprefix = outprefix,
    attr_ID_field = "GOID"
  )
  return(filtered_result)
}


# Filter annotation table by minimum number of genes.
# Input:  gene_geneset_table: annotation table in the format col:(SYMBOL, attrb1, attrb2) row: the membership of each gene . min_ngene: min number of gene
# max_ngene: max number of gene to be included in the output, geneset_description: dataframe of description of geneset)
filter_annotation_by_ngene = function(gene_geneset_table,
                                      geneset_description,
                                      min_ngene = 5,
                                      max_ngene = 2000,
                                      outprefix = paste0(data_path, "/GO"),
                                      attr_ID_field = "GOID")
{
  gene_count = data.frame(n_gene = colSums(gene_geneset_table[, c(2:ncol(gene_geneset_table))]))
  colnames(gene_geneset_table)[1] = "SYMBOL"
  gene_count$attr_name = rownames(gene_count)
  selected_attrs = gene_count[(gene_count$n_gene >= min_ngene &
                                 gene_count$n_gene <= max_ngene), "attr_name"]
  gene_geneset_table = gene_geneset_table[, c("SYMBOL", selected_attrs)]
  geneset_description = geneset_description[geneset_description[, attr_ID_field] %fin% selected_attrs,]
  write.table(
    gene_geneset_table,
    file = paste(outprefix, "_gene_geneset_table.csv", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    geneset_description,
    file = paste(outprefix, "_geneset_description.csv", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  return(list(gene_geneset_table, geneset_description))
}












# Select top n_token frequent words  of each annotation group in annotation_group_field for final annoatation
get_word_annotation = function(subnet_annotation_table = sn_an,
                               gene_geneset_table,
                               annotation_def = MSigDB_description,
                               definition_field = "DESCRIPTION_BRIEF",
                               annotation_id_field = "SYSTEMATIC_NAME",
                               n_token = 10,
                               annotation_group_field = "ONTOLOGY")
{
  annotation_terms = colnames(subnet_annotation_table)[2:ncol(subnet_annotation_table)]
  annotation_def = annotation_def[annotation_def[, annotation_id_field] %in% annotation_terms, ] # filter definition table
  
  ngene_per_geneset = data.frame(n_genes = colSums(gene_geneset_table[2:ncol(gene_geneset_table)]),
                                 check.names = FALSE)
  ngene_per_geneset$annotation_id = rownames(ngene_per_geneset)
  annotation_def = merge(
    annotation_def,
    ngene_per_geneset,
    by.x = annotation_id_field,
    by.y = "annotation_id",
    all.x = TRUE
  )
  if (!is.null(annotation_group_field))
  {
    annotation_terms = colnames(subnet_annotation_table)[2:ncol(subnet_annotation_table)]
    annotation_groups = unique(annotation_def[, annotation_group_field])
  } else
    annotation_groups = c("all")
  subnet_names = as.character(subnet_annotation_table$subnet_name)
  subnet_word_annotation = data.frame(subnet_name = subnet_names)
  for (group in  annotation_groups)
  {
    #subnet_word_annotation$annotation = ""
    # for (subnet_i in c(2:length(subnet_names)))  # skipping first row as this is the root subnet
    if (!is.null(annotation_group_field))
    {
      group_terms = annotation_def[annotation_def[, annotation_group_field] == group, annotation_id_field]
    } else
    {
      group_terms = annotation_def[, annotation_id_field]
    }
    for (subnet_i in c(1:length(subnet_names)))
    {
      #annotated_pathways =subnet_annotation_table[subnet_i, c(4:length(subnet_annotation_table))]
      annotated_pathways = subnet_annotation_table[subnet_i, group_terms]
      annotated_pathways = t(annotated_pathways)
      annotated_pathways = data.frame(annotated_pathways)
      colnames(annotated_pathways)[1] = "path_activation"
      annotated_pathways$pathway = rownames(annotated_pathways)
      annotated_pathways = data.frame(annotated_pathways[annotated_pathways$path_activation ==
                                                           1, ])
      annotated_pathways = merge(annotated_pathways,
                                 annotation_def,
                                 by.x = "pathway",
                                 by.y = annotation_id_field)
      word_annotation = annotated_pathways[, definition_field]
      token_annotation = select_top_words(word_annotation, excluded_word_list, n_token)
      
      subnet_word_annotation[subnet_i, group] = token_annotation
      subnet_name = subnet_names[subnet_i]
      #print(paste(group, subnet_name, token_annotation))
    }
    #print(paste("finish this group", group))
  }
  
  return(subnet_word_annotation)
}

# from the vector of annotation , select top most frequently used words
select_top_words = function(word_annotation,
                            excluded_word_list,
                            n_token = 10)
{
  word_annotation = as.vector(word_annotation)
  word_annotation = strsplit(word_annotation, " ")
  tokens = unlist(word_annotation)
  tokens = data.frame(token = tokens)
  
  if (nrow(tokens) > 1)
  {
    token_count = ddply(tokens, ~ token, nrow)
    token_count = token_count[!(token_count$token %in% excluded_word_list),]
    
    token_count = token_count[order(token_count$V1, decreasing = TRUE), ]
    if (nrow(token_count) > n_token)
    {
      token_count = as.character(token_count[1:n_token, "token"])
    } else
    {
      token_count = as.character(token_count[, "token"])
    }
    
    token_annotation = paste(token_count, collapse = ", ")
    
  }
  else
  {
    token_annotation = "-"
  }
  return(token_annotation)
  
}




#parse MSigDB xml file
#min_ngene: mimimul number of gene in pathway to be parsed
parse_MSigDB = function(xml_file = MSigDB_file,
                        outdir = paste0(data_path),
                        min_ngene = 5,
                        max_ngene = 2000,
                        organism = "Homo sapiens",
                        excluding_GO = TRUE,
                        canonical_pathway_only = TRUE)
{
  require(XML)
  input_data <- xmlParse(xml_file)
  xml_data <- xmlToList(input_data)
  gs_description = data.frame()
  for (i in c(1:length(xml_data)))
  {
    #print(i)
    record = xml_data[[i]]
    fields = names(record)
    for (field in fields)
    {
      gs_description[i, field] = record[field]
    }
    
  }
  gs_description = gs_description[gs_description$ORGANISM == organism,]
  if (excluding_GO)
  {
    gs_description = gs_description[gs_description$CONTRIBUTOR_ORG != "GO", ]
  }
  
  if (canonical_pathway_only)
  {
    gs_description = gs_description[(
      gs_description$SUB_CATEGORY_CODE %in% c("CP:BIOCARTA", "CP:KEGG", "CP:REACTOME")
    ) |
      (gs_description$CATEGORY_CODE == "H")  &
      !is.na(gs_description$STANDARD_NAME), ]
  }
  
  genesets = gs_description$SYSTEMATIC_NAME
  
  n_geneset = nrow(gs_description)
  gene_list = list()
  geneset_list =  list()
  for (i in c(1:n_geneset))
  {
    geneset = genesets[i]
    print(paste("parsing pathway", i, geneset))
    genes = gs_description[i, "MEMBERS_SYMBOLIZED"]
    genes = unlist(strsplit(genes, ","))
    n_gene = length(genes)
    if (n_gene >= min_ngene & n_gene < max_ngene)
    {
      geneset_list[[i]] = geneset  # only select geneset that satisfies requirement of min, max n of genes
      gene_list[[i]] = genes
    }
  }
  gs_description  = gs_description[gs_description$SYSTEMATIC_NAME %in%  geneset_list,]
  unique_genes = unique(unlist(gene_list))
  gene_geneset_df = data.frame(SYMBOL = unique_genes)
  
  for (i in c(1:length(geneset_list)))
  {
    print(paste("inserting pathway", i, geneset))
    geneset = geneset_list[[i]]
    genes = gene_list[[i]]
    require(fastmatch)
    genes_match = fmatch(genes, unique_genes)
    gene_geneset_df[, geneset] = 0
    gene_geneset_df[genes_match, geneset] = 1
  }
  
  geneset_summary = unlist(colSums(gene_geneset_df[2:ncol(gene_geneset_df)]))
  gs_description$n_gene = geneset_summary
  write.table(
    gene_geneset_df,
    file = paste(outdir, "/MSigDB_gene_geneset_table.csv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    gs_description,
    file = paste(outdir, "/MSigDB_geneset_description.csv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  return(list(gene_geneset_df, df))
}

# get list of annotation table (annotation_term, subnet1, subnet2, subnet, ...) and summary table (subnet_name, n_annotation)
geneset_cluster_annotate = function(subnet_annotation_table,
                                    annotation_def,
                                    definition_field = "DESCRIPTION_BRIEF",
                                    annotation_id_field = "SYSTEMATIC_NAME",
                                    gene_geneset_table,
                                    subnet_qvalues,
                                    subnets)
{
  sn_an_final = data.frame(t(subnet_annotation_table[, c(2:ncol(subnet_annotation_table))]), check.names = FALSE)
  sn_an_sum = data.frame(n_annotations = colSums(sn_an_final),
                         check.names = FALSE)
  sn_an_sum$subnet = rownames(sn_an_sum)
  sn_an_final = cbind(data.frame(
    annotation_id = rownames(sn_an_final),
    check.names = FALSE
  ),
  sn_an_final)
  
  sn_an_by_coverage = sn_an_final
  subnet_names = colnames(sn_an_by_coverage)
  sn_an_by_coverage[, "annotation_id"] = as.character(sn_an_by_coverage[, "annotation_id"])
  
  for (i in c(1:nrow(sn_an_by_coverage)))
    # loop for each annotation term
  {
    geneset_name  = sn_an_by_coverage[i, "annotation_id"]
    geneset_genes = gene_geneset_table[gene_geneset_table[, geneset_name] == 1, 1]
    #print(geneset_name)
    for (j in c(2:ncol(sn_an_by_coverage)))
      # loop for each subnet
    {
      #is_active = sn_an_by_coverage[i, j]
      is_active = .subset2(sn_an_by_coverage, j)[i]
      if (is_active == 1)
        # if the term is active for this subnet, then calculate the number of overlapping genes between the term and the subnet
      {
        subnet = subnets[[j]]
        subnet_name = graph_attr(subnet, name = "graph_name")
        subnet_genes = graph_attr(subnet, name = "element_node_names")
        #subnet_genes_match = fmatch(subnet_genes, geneset_genes$SYMBOL)
        #subnet_at = at[subnet_genes_match,2]
        share_genes = geneset_genes[geneset_genes %fin% subnet_genes]
        coverage = length(share_genes) / length(subnet_genes)
        sn_an_by_coverage[i, j] = coverage
      }
    }
    
  }
  
  
  ngene_per_geneset = data.frame(n_genes = colSums(gene_geneset_table[2:ncol(gene_geneset_table)]),
                                 check.names = FALSE)
  
  ngene_per_geneset$annotation_id = rownames(ngene_per_geneset)
  annotation_def = merge(
    annotation_def,
    ngene_per_geneset,
    by.x = annotation_id_field,
    by.y = "annotation_id",
    all.x = TRUE
  )
  sn_an_final = merge(
    annotation_def,
    sn_an_final,
    by.x = annotation_id_field,
    by.y = "annotation_id",
    all.y = TRUE
  )
  # Calculate the number of gene that share between a subnet and a term.
  
  return(list(sn_an_final, sn_an_sum, sn_an_by_coverage))
}

get_one_cluster_annotation = function(cluster_annotations,
                                      cluster_name = "1.7",
                                      selected_field = c("SYSTEMATIC_NAME", "STANDARD_NAME", "DESCRIPTION_BRIEF"))
{
  #cluster_name = paste("X",cluster_name,sep="")
  pathway = cluster_annotations[, c(selected_field, cluster_name)]
  selected_pathway = pathway[pathway[, cluster_name] == 1, ]
  return(selected_pathway)
}


# get gene summary , input is Entrez ID. need to convert symbol to Entrez ID

get_gene_summary = function(gene_ids = c(14912, 712))
{
  library(rentrez)
  genes <- entrez_summary(db = "gene", gene_ids)
  out = data.frame(ENTREZID = gene_ids)
  out$description = ""
  out$summary = ""
  for (i in c(1:length(gene_ids)))
  {
    description = genes[[i]]$description
    summary = genes[[i]]$summary
    out[i, "summary"] = summary
    out[i, "description"] = description
  }
  return(out)
}

get_all_gene_info = function(organism = "Homo sapiens",
                             batch_size = 300,
                             outfile = paste0(data_path, "/all_gene_summary.csv"))
{
  library(AnnotationHub)
  ah = AnnotationHub()
  orgs = subset(ah, ah$rdataclass == "OrgDb")
  org = query(orgs, organism)[[1]]
  all_genes =  unique(keys(org, keytype = "ENTREZID"))
  genes_info = select(
    org,
    keys = all_genes,
    columns = c("ENTREZID", "SYMBOL", "PFAM", "UNIPROT"),
    keytype = "ENTREZID"
  )
  first_gene_Id = 1
  total_genes = length(all_genes)
  
  while (first_gene_Id <= total_genes)
  {
    last_gene_Id = first_gene_Id + batch_size - 1
    print(paste("getting gene info from", first_gene_Id, "to", last_gene_Id))
    if (last_gene_Id > total_genes)
    {
      last_gene_Id = total_genes
    }
    gene_summary = get_gene_summary(all_genes[first_gene_Id:last_gene_Id])
    if (exists("combined_gene_summary"))
    {
      combined_gene_summary = rbind(combined_gene_summary, gene_summary)
    } else
    {
      combined_gene_summary = gene_summary
    }
    first_gene_Id = last_gene_Id + 1
  }
  
  combined_gene_summary = merge(
    genes_info,
    combined_gene_summary,
    by.x = "ENTREZID",
    by.y = "ENTREZID",
    all.x = TRUE,
    sort = FALSE
  )
  write.table(
    combined_gene_summary,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  return(combined_gene_summary)
}

cluster_annotate_by_gene_description = function(subnets,
                                                gene_description_table,
                                                gene_id_field = "SYMBOL",
                                                description_field = "summary",
                                                n_token = 10,
                                                n_top_frequent_word = 200)
{
  all_word_annotation = gene_description_table[, description_field]
  top_repeated_tokens = select_top_words(
    word_annotation = all_word_annotation,
    excluded_word_list = c(),
    n_token = n_top_frequent_word
  )
  top_repeated_tokens = strsplit(top_repeated_tokens, split = ", ")
  excluded_word_list = c(top_repeated_tokens, excluded_word_list)
  for (subnet in subnets)
  {
    subnet_name = graph_attr(subnet, name = "graph_name")
    genes = graph_attr(subnet, name = "element_node_names")
    gene_descriptions = gene_description_table[gene_description_table[, gene_id_field] %in% genes, c(gene_id_field, description_field)]
    gene_descriptions = unique(gene_descriptions)
    word_annotation = gene_descriptions[, description_field]
    
    word_annotation = select_top_words(
      word_annotation = word_annotation,
      excluded_word_list = excluded_word_list,
      n_token = n_token
    )
    #print(paste(subnet_name, word_annotation))
    subnet_annotation = data.frame(subnet_name = subnet_name, annotation = word_annotation)
    if (exists("combined_subnet_annotation"))
    {
      combined_subnet_annotation = rbind(combined_subnet_annotation, subnet_annotation)
    } else
    {
      combined_subnet_annotation = subnet_annotation
    }
  }
  
  return(combined_subnet_annotation)
}

## annotation sources is a list of annotation_source. Each annotation_source is a list(gene_geneset_table , geneset_description, definition_field, annotation_id_field )
get_annotation_sources = function(excluded_terms = NULL)
{
  if (!file.exists(paste0(data_path, "/GO_gene_geneset_table.csv")))
  {
    hugo =  get_human_go()
    rm(hugo)
  }
  if (!file.exists(paste0(data_path, "/MSigDB_gene_geneset_table.csv")))
  {
    MSigDB = parse_MSigDB()  # generate MSingDB gene_geneset_table and description table if not exists
    rm(MSigDB)
  }
  MSigDB_df <-
    read.delim(
      paste0(data_path, "/MSigDB_gene_geneset_table.csv"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  MSigDB_description <-
    read.delim(
      paste0(data_path, "/MSigDB_geneset_description.csv"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  go_gene_geneset_table = read.delim(
    paste0(data_path, "/GO_gene_geneset_table.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  go_geneset_description = read.delim(
    paste0(data_path, "/GO_geneset_description.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(excluded_terms))
  {
    go_geneset_description = go_geneset_description[!go_geneset_description$TERM %in% excluded_GO_terms, ]
    selected_goid = go_geneset_description$GOID
    go_gene_geneset_table = go_gene_geneset_table[, c("SYMBOL", selected_goid)]
  }
  
  go = list(
    gene_geneset_table = go_gene_geneset_table,
    geneset_description = go_geneset_description,
    definition_field = "TERM",
    annotation_id_field = "GOID",
    annotation_group_field = "ONTOLOGY"
  )
  MSigDB = list(
    gene_geneset_table = MSigDB_df,
    geneset_description = MSigDB_description,
    definition_field = "DESCRIPTION_BRIEF",
    annotation_id_field = "SYSTEMATIC_NAME",
    annotation_group_field = NULL
  )
  annotation_sources = list(go, MSigDB)
  return(annotation_sources)
}

# Combine GO, MSigDB, gene definition annotation
combined_word_annotations = function(annotation_score_list,
                                     annotation_sources,
                                     subnets,
                                     gene_description_file = paste0(data_path, "/all_gene_summary.csv"),
                                     n_top_annotation = 1,
                                     importance_field = "score",
                                     n_top_genes = 3,
                                     graph)
{
  print(paste("n_top_genes", n_top_genes))
  combined_annotations = list()
  for (i in c(1:length(annotation_sources)))
  {
    node_attr_score_df = annotation_score_list[[i]]
    out  = get_subnet_enrichment_score_from_subnets_list(
      subnets = subnets,
      node_attrs = node_attr_score_df,
      sampling_size = 1000,
      use_random_sampling = FALSE
    )
    subnet_pvalues = out[[1]]
    subnet_qvalues =  out[[2]]
    subnet_enrichment_values = out[[3]]
    subnet_qvalues[is.na(subnet_qvalues)] = 1
    annotation_source = annotation_sources[[i]]
    sn_an = subnet_annotation_level_by_level(
      subnet_qvalues = subnet_qvalues,
      subnets = subnets,
      gene_geneset_table = annotation_source$gene_geneset_table
    )
    
    #convert from annotation id to word annotation for each subnet
    cluster_annotations = geneset_cluster_annotate(
      subnet_annotation_table = sn_an,
      annotation_def = annotation_source$geneset_description,
      definition_field = annotation_source$definition_field,
      annotation_id_field =
        annotation_source$annotation_id_field ,
      gene_geneset_table = annotation_source$gene_geneset_table,
      subnet_qvalues = subnet_qvalues,
      subnets = subnets
    )
    cl_annotation = cluster_annotations[[1]]
    cl_annotation_sum = cluster_annotations[[2]]
    if (!exists("cl_annotation_by_coverage"))
    {
      cl_annotation_by_coverage = cluster_annotations[[3]]
    } else
    {
      cl_annotation_by_coverage = rbind(cl_annotation_by_coverage, cluster_annotations[[3]])
    }
    combined_annotations[[i]] = cluster_annotations
    if (!exists("combined_subnet_qvalues"))
    {
      combined_subnet_qvalues = subnet_qvalues
    } else
    {
      combined_subnet_qvalues = fast_merge(
        combined_subnet_qvalues,
        subnet_qvalues[, c(1, 4:ncol(subnet_qvalues))],
        by.x = "subnet_name",
        by.y = "subnet_name",
        all.x = TRUE
      )
    }
    
  }
  
  # get combined term definition table
  
  for (as_index in c(1:length(annotation_sources)))
  {
    #print(as_index)
    annotation_source = annotation_sources[[as_index]]
    geneset_description = annotation_source$geneset_description
    definition_field = annotation_source$definition_field
    annotation_id_field = annotation_source$annotation_id_field
    id_definition_table = geneset_description[, c(annotation_id_field, definition_field)]
    colnames(id_definition_table) = c("annotation_id", "annotation")
    if (exists("combined_id_definition_table"))
    {
      combined_id_definition_table = rbind(combined_id_definition_table, id_definition_table)
    } else
    {
      combined_id_definition_table = id_definition_table
    }
  }
  
  # Select top annotation for each subnet
  cl_top_annotation = data.frame(subnet_name = colnames(cl_annotation_by_coverage)[2:ncol(cl_annotation_by_coverage)])
  for (j in c(2:ncol(cl_annotation_by_coverage)))
  {
    if (j == ncol(cl_annotation_by_coverage))
    {
      print("stop selecting top annotation for each subnet")
    }
    one_cluster = cl_annotation_by_coverage[, c(1, j)]
    one_cluster = one_cluster[one_cluster[, 2] > 0, ]
    if (nrow(one_cluster) > 0)
    {
      one_cluster = one_cluster[order(one_cluster[, 2], decreasing = TRUE), ]
      
      if (nrow(one_cluster) > n_top_annotation)
      {
        top_annotations = one_cluster[1:n_top_annotation, 1]
      } else
      {
        top_annotations = one_cluster[, 1]
      }
      top_annotations_string = paste(top_annotations, collapse = "\n")
      top_annotations_df = data.frame(annotation_id = top_annotations)
      top_annotations_df = fast_merge(
        top_annotations_df,
        combined_id_definition_table,
        by.x = "annotation_id",
        by.y = "annotation_id",
        all.x = TRUE,
        all.y = FALSE
      )
      combined_top_def = paste(top_annotations_df$annotation, collapse = "\n")
      #print(combined_top_def)
    } else
    {
      combined_top_def = ""
      top_annotations_string = ""
    }
    cl_top_annotation[j - 1 , "top_annotation"] = top_annotations_string
    cl_top_annotation[j - 1 , "annotation"] = combined_top_def
    
  }
  
  gene_attr = data.frame(vertex_attr(graph))
  for (i in c(2:length(subnets)))
  {
    level = graph_attr(subnets[[i]], name = "level")
    cl_top_annotation[i - 1 , "level"] = level
    
    
    subnet = subnets[[i]]
    subnet_genes = graph_attr(subnet, name = "element_node_names")
    get_top_genes = function(field)
    {
      subnet_genes_most_important = gene_attr[gene_attr$name %fin% subnet_genes, c("name", field)]
      #print(paste(2244, subnet_genes_most_important))
      subnet_genes_most_important = subnet_genes_most_important[subnet_genes_most_important[, field] > 0,]
      #print(paste(2246, subnet_genes_most_important))
      subnet_genes_most_important = subnet_genes_most_important[!is.na(subnet_genes_most_important$name),]
      #print(paste(2248, subnet_genes_most_important))
      subnet_genes_most_important = subnet_genes_most_important[order(subnet_genes_most_important[, field], decreasing = TRUE), ]
      #print(2250)
      #print(subnet_genes_most_important)
      #print(nrow(subnet_genes_most_important))
      #print(n_top_genes)
      #print(2254)
      if (nrow(subnet_genes_most_important) > n_top_genes)
      {
        subnet_genes_most_important = subnet_genes_most_important[1:n_top_genes, "name"]
      } else
      {
        subnet_genes_most_important = subnet_genes_most_important[, "name"]
      }
      subnet_genes_most_important = paste(subnet_genes_most_important, collapse = ", ")
      return(subnet_genes_most_important)
    }
    cl_top_annotation[i - 1 , "top_hubs"] =  get_top_genes("degree")
    if (!is.null(importance_field))
    {
      cl_top_annotation[i - 1 , "top_master_regulators"] =  get_top_genes(importance_field)
    }
    
  }
  
  
  
  
  # Add annotation column to cl_annotation_by_coverage
  
  cl_annotation_by_coverage = fast_merge(
    cl_annotation_by_coverage,
    combined_id_definition_table,
    by.x = "annotation_id",
    by.y = "annotation_id",
    all.x = TRUE
  )
  cl_annotation_by_coverage = cl_annotation_by_coverage[, c(1, ncol(cl_annotation_by_coverage), c(2:(ncol(
    cl_annotation_by_coverage
  ) - 1)))]
  return(
    list(
      combined_annotations = combined_annotations,
      cl_top_annotation = cl_top_annotation,
      cl_annotation_by_coverage = cl_annotation_by_coverage,
      combined_subnet_qvalues = combined_subnet_qvalues
    )
  )
}



test_multiple_annotation_datasets = function()
{
  network_file = "../data/nha_gsc_net_nbhscore.graphml"
  graph_format = "graphml"
  #outfile = "nha_gsc_net_openord3D_2D.csv"
  clustering_method = "louvain"
  max_cluster_size = 100
  options(expressions = 10000)
  ndim = 3
  layout_method = "fr"
  space_ratio = 2
  node_score_field = "logFC"
  enrichment_score_field = c("node_smooth_score")
  sampling_size = 1000
  beta = 0.5
  n_top_annotation = 3
  importance_field = "score"
  n_top_genes = 3
  
  graph = read_graph(file = network_file, format = graph_format)
  graph = as.undirected(graph, mode = c("collapse"))
  annotation_sources = get_annotation_sources(excluded_terms = excluded_GO_terms)
  attr_list = list()
  for (i in c(1:length(annotation_sources)))
  {
    attr_list[[i]] = annotation_sources[[i]]$gene_geneset_table
  }
  
  
  result = get_network_layout(
    graph = graph,
    clustering_method = "louvain",
    max_cluster_size = 50,
    ndim = 3,
    layout_method = "fr",
    space_ratio = 2,
    node_score_field = "logFC",
    beta = 0.5,
    attr = attr_list
  )
  
  
  layout = result[[1]]
  subnets = result[[2]]
  annotation_score_list = result[[3]]
  
  combined_annotations = combined_word_annotations(
    annotation_score_list,
    annotation_sources,
    subnets = subnets,
    n_top_annotation = 1
  )
  cl_top_annotation = combined_annotations[[2]]
  
  
  return(cl_top_annotation)
}





# Extract subnets at all levels that the gene belong to and all the annotation associated with these subnets
#Input: gene_name: gene name; layout: final layout, result from  get_network_layout_pipeline. subnets_attr: description of each subnet,  result from get_network_layout_pipeline
extract_gene_function = function(gene_name = "SNAI1",
                                 layout,
                                 subnets_attr,
                                 outdir = "../data/",
                                 subnet_annotation_by_coverage,
                                 combined_subnet_qvalues,
                                 combined_gene_geneset_table,
                                 diff_exp_field = "logFC",
                                 selected_fields = c(
                                   "node_name",
                                   "score",
                                   diff_exp_field,
                                   "description",
                                   "summary",
                                   "subnet_name",
                                   "degree",
                                   "pval",
                                   "fdr",
                                   "annotation"
                                 ),
                                 draw_heatmap = TRUE,
                                 gene_expression_table = NULL,
                                 sample_group_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_sample_description.csv",
                                 displayed_groups = NULL,
                                 n_top_genes = 50,
                                 importance_field = "score",
                                 is_log_gene_expression_table = FALSE,
                                 n_top_pathway_to_draw_heatmap = 5,
                                 baseline_samples = NULL,
                                 annotation_sources = NULL,
                                 use_cluster = FALSE,
                                 baseline_table = NULL)
{
  cut_off_field = diff_exp_field
  subnet_name = as.character(layout[layout$node_name == gene_name, "subnet_name"])
  # Getting all parent subnets for the gene
  subnet_name_parsed = unlist(strsplit(subnet_name, split = ".", fixed = TRUE))
  parent_subnet_name =  "1"
  n_level = length(subnet_name_parsed)
  subnet_level_names = c()
  for (i in c(2:n_level))
  {
    subnet_level_name = paste(parent_subnet_name, subnet_name_parsed[i], sep =
                                ".")
    subnet_level_names = c(subnet_level_names, subnet_level_name)
    parent_subnet_name = subnet_level_name
  }
  gene_subnets = subnets_attr[subnets_attr$graph_name %in% subnet_level_names, ]
  for (i in c(1:nrow(gene_subnets)))
  {
    split_top_annotation = unlist(strsplit(gene_subnets[i, "top_annotation"], split =
                                             "\n"))
    gene_subnets[i, "top_annotation"] = paste(split_top_annotation, collapse = "; ")
    split_annotation = unlist(strsplit(gene_subnets[i, "annotation"], split =
                                         "\n"))
    gene_subnets[i, "annotation"] = paste(split_annotation, collapse = "; ")
  }
  if (!is.null(outdir))
  {
    dir.create(outdir)
    outdir = paste(outdir, "/", gene_name, sep = "")
    dir.create(outdir)
    outfile = paste(outdir, "/", gene_name, "_function.txt", sep = "")
    write.table(
      gene_subnets,
      file = outfile,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    for (i in c(1:nrow(gene_subnets)))
    {
      subnet = as.character(gene_subnets$graph_name[i])
      subnet_data = extract_subnet_data(
        final_layout = layout,
        subnet = subnet,
        subnets_attr = subnets_attr,
        subnet_annotation_by_coverage = subnet_annotation_by_coverage,
        combined_subnet_qvalues = combined_subnet_qvalues,
        combined_gene_geneset_table = combined_gene_geneset_table,
        gene_description_file = paste0(data_path, "/all_gene_summary.csv"),
        qvalue_cutoff = 0.05,
        outdir = paste0(outdir,"/",subnet ),
        cut_off_field = cut_off_field ,
        sort_fields = c(cut_off_field , "score"),
        cut_off_value = 1,
        selected_fields = selected_fields,
        draw_heatmap = draw_heatmap,
        gene_expression_table = gene_expression_table,
        sample_group_table = sample_group_table,
        displayed_groups = displayed_groups,
        n_top_genes = n_top_genes,
        importance_field = importance_field,
        is_log_gene_expression_table = is_log_gene_expression_table,
        n_top_pathway_to_draw_heatmap =  n_top_pathway_to_draw_heatmap ,
        baseline_samples = baseline_samples,
        annotation_sources = annotation_sources,
        use_cluster = use_cluster,
        baseline_table = baseline_table
      )
      
    }
    
    
    
  }
  return(gene_subnets)
}







extract_gene_function_from_network_layout_pipeline_after_saving_RDS = function(gene_name = "SNAI1",
                                                                               network_layout_pipeline_after_saving_RDS = "get_network_layout_pipeline_after_saving.RDS",
                                                                               layout,
                                                                               subnets_attr,
                                                                               outdir = "../data/",
                                                                               subnet_annotation_by_coverage,
                                                                               combined_subnet_qvalues,
                                                                               combined_gene_geneset_table,
                                                                               diff_exp_field = "logFC",
                                                                               selected_fields = c(
                                                                                 "node_name",
                                                                                 "score",
                                                                                 diff_exp_field,
                                                                                 "description",
                                                                                 "summary",
                                                                                 "subnet_name",
                                                                                 "degree",
                                                                                 "pval",
                                                                                 "fdr",
                                                                                 "annotation"
                                                                               ),
                                                                               baseline_samples = NULL)
{
  #' extract gene function from network layout pipeline data after saving
  #'
  #' Using saved RDS to extract gene function. RDS is the output of save_network_layout function
  #'
  layout = network_layout_pipeline_after_saving_RDS$layout
  subnets_attr = network_layout_pipeline_after_saving_RDS$subnets_attr
  subnet_annotation_by_coverage =  network_layout_pipeline_after_saving_RDS$subnet_annotation_by_coverage
  combined_subnet_qvalues = network_layout_pipeline_after_saving_RDS$combined_subnet_qvalues
  combined_gene_geneset_table = network_layout_pipeline_after_saving_RDS$combined_gene_geneset_table
  
  
  out = extract_gene_function_from_network_layout_pipeline_after_savingRDS(
    gene_name = gene_name,
    layout = layout,
    subnets_attr = subnets_attr,
    outdir = outdir,
    subnet_annotation_by_coverage =  subnet_annotation_by_coverage,
    combined_subnet_qvalues = combined_subnet_qvalues,
    combined_gene_geneset_table = combined_gene_geneset_table,
    diff_exp_field = diff_exp_field,
    selected_fields = selected_fields,
    baseline_samples = baseline_samples
  )
  return(out)
  
}



get_hierarchy_network = function(final_layout,
                                 subnets_attr,
                                 outdir = "../data/layout",
                                 diff_exp_field = "logFC")
{
  gene_edges = data.frame(node_1 = final_layout$subnet_name,
                          node_2 = final_layout$node_name)
  subnet_names = as.character(subnets_attr$graph_name[2:nrow(subnets_attr)])
  subnet_parent_names = c()
  for (subnet_name in subnet_names)
  {
    parent_name = get_subnet_parent_name(subnet_name)
    subnet_parent_names = c(subnet_parent_names, parent_name)
  }
  subnet_edges = data.frame(node_1 = subnet_parent_names, node_2 = subnet_names)
  combined_edges = rbind(subnet_edges, gene_edges)
  #print(diff_exp_field)
  selected_fields = c(
    "node_name",
    "level",
    "net_type",
    "subnet_score",
    "degree",
    diff_exp_field ,
    "subnet_qvalue",
    "enrichment_subnet_score_with_sign"
  )
  #print(2392)
  #print(selected_fields)
  #print(colnames(final_layout))
  genes_attr_short = final_layout[, selected_fields]
  genes_attr_short$n_genes = 1
  subnets_attr_short = subnets_attr[, c(
    "graph_name",
    "level",
    "net_type",
    "n_genes",
    "subnet_score",
    "subnet_qvalue",
    "enrichment_subnet_score_with_sign",
    "annotation",
    "top_hubs",
    "top_master_regulators"
  )]
  colnames(subnets_attr_short)[1] = "node_name"
  hierarchy_network_attrs = rbind.fill(subnets_attr_short, genes_attr_short)
  hierarchy_network_attrs[, "combined_node_name"] = paste(
    hierarchy_network_attrs$node_name,
    "\n",
    hierarchy_network_attrs$top_master_regulators,
    "\n",
    hierarchy_network_attrs$annotation
  )
  hierarchy_network_attrs[hierarchy_network_attrs$net_type == 3, "combined_node_name"] = as.character(hierarchy_network_attrs[hierarchy_network_attrs$net_type ==
                                                                                                                                3, "node_name"])
  hierarchy_network_attrs[hierarchy_network_attrs$level == 1, "combined_node_name"] = "1 \n ROOT"
  hierarchy_network = graph_from_data_frame(combined_edges, directed = TRUE, vertices = hierarchy_network_attrs)
  dir.create(outdir)
  graph_file = paste(outdir, "/hierarchy_graph.graphml", sep = "")
  #layout = layout_with_lgl(hierarchy_network)
  #plot(hierarchy_network, layout = layout, vertex.size =0.1, vertex.label.cex =0.1, edge.arrow.size=0.01)
  write_graph(hierarchy_network, file = graph_file, format = "graphml")
  
}

get_subnet_parent_name = function(subnet_name)
{
  dot_pos = as.vector(gregexpr(".", subnet_name, fixed = TRUE)[[1]])
  if (length(dot_pos) > 0)
  {
    last_pos = dot_pos[length(dot_pos)]
    parent_name = substring(subnet_name, 1, last_pos - 1)
  } else
  {
    parent_name = ""
  }
  return(parent_name)
}

# add external compartment to final layout data frame
get_external_compartment_field = function(final_layout,
                                          GO_full_file =
                                            paste0(data_path, "/GO_full_gene_geneset_table.csv"))
{
  if (file.exists(GO_full_file))
  {
    go_gene_geneset_table = read.delim(GO_full_file,
                                       stringsAsFactors = FALSE,
                                       check.names = FALSE)
    go_geneset_description = read.delim(
      paste0(data_path, "/GO_full_geneset_description.csv"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else
  {
    go = get_human_go(
      min_ngene = 1,
      max_ngene = 500000,
      outprefix = paste0(data_path, "/GO_full")
    )
    go_gene_geneset_table = go[[1]]
    go_geneset_description = go[[2]]
  }
  
  go_geneset_description_CC =  go_geneset_description[go_geneset_description$ONTOLOGY ==
                                                        "CC", ]
  go_gene_geneset_table_CC = go_gene_geneset_table[, c("SYMBOL", go_geneset_description_CC$GOID)]
  external_compartment = c(
    "external side of plasma membrane",
    "basolateral plasma membrane",
    "apical plasma membrane",
    "plasma membrane",
    "lateral plasma membrane",
    "apicolateral plasma membrane",
    "extrinsic component of plasma membrane",
    "extrinsic component of membrane",
    "plasma membrane protein complex",
    "plasma membrane",
    "cell surface",
    "extracellular region",
    "extracellular region part",
    "extracellular exosome complex",
    "extracellular matrix",
    "extracellular space"
  )
  GO_external_ID = go_geneset_description_CC[go_geneset_description_CC$TERM %in% external_compartment,]
  go_gene_geneset_table_external =   go_gene_geneset_table_CC[, c("SYMBOL", GO_external_ID$GOID)]
  go_gene_geneset_table_external$external_sum = apply(go_gene_geneset_table_external[2:ncol(go_gene_geneset_table_external)], 1, sum)
  go_gene_geneset_table_external = go_gene_geneset_table_external[go_gene_geneset_table_external$external_sum > 0,]
  external_genes =  as.character(go_gene_geneset_table_external$SYMBOL)
  final_layout$is_external = 0
  final_layout[final_layout$node_name %in% external_genes, "is_external"] =
    1
  return(final_layout)
}

external_compartment_analysis = function(final_layout,
                                         cut_off_field = "logFC",
                                         cut_off_value = 1)
{
  if (is.null(final_layout$is_external[1]))
  {
    final_layout = get_external_compartment_field(final_layout)
  }
  genes_external_signif_increase = final_layout[final_layout$is_external == 1 &
                                                  final_layout$subnet_qvalue < 0.05 &
                                                  final_layout[, cut_off_field] > cut_off_value, ]
  genes_external_signif_increase =  genes_external_signif_increase[order(genes_external_signif_increase[, cut_off_field], decreasing = TRUE), ]
  genes_external_signif_increase_annotated_short = genes_external_signif_increase[, c("node_name",
                                                                                      "subnet_name",
                                                                                      cut_off_field,
                                                                                      "description",
                                                                                      "summary",
                                                                                      "annotation")]
  return(
    list(
      genes_external_signif_increase,
      genes_external_signif_increase_annotated_short
    )
  )
}




extract_significant_genes_of_a_cluster = function(final_layout,
                                                  cluster = "1.5.6",
                                                  diff_cutoff = 1,
                                                  qvalue_cutoff = 0.05)
{
  sign_genes = final_layout[final_layout$subnet_qvalue < qvalue_cutoff &
                              substr(final_layout$subnet_name, 1, nchar(cluster)) == cluster &
                              abs(final_layout$diff) > diff_cutoff,]
  return(sign_genes)
}





# Create heatmap of genes belonging to subnet and its children subnetworks
# gene_expression_table: first column: gene name, other columns: sample ID. Row is the gene expression value.

subnet_heatmap = function(final_layout,
                          subnet = "1.5.6",
                          diff_cutoff = 1,
                          qvalue_cutoff = 0.05,
                          gene_expression_table)
{
  
}


# Create subnet network for export in Paraview or Cytoscape with node is a gene.

extract_subnet_data = function(final_layout,
                               subnet = "1.5.5.3",
                               subnets_attr,
                               subnet_annotation_by_coverage,
                               combined_subnet_qvalues,
                               combined_gene_geneset_table,
                               gene_description_file = paste0(data_path, "/all_gene_summary.csv"),
                               qvalue_cutoff = 0.05,
                               outdir = "../data/layout",
                               cut_off_field = "logFC" ,
                               sort_fields = c(cut_off_field , "score"),
                               cut_off_value = 1,
                               selected_fields = c(
                                 "node_name",
                                 "score",
                                 cut_off_field ,
                                 "description",
                                 "summary",
                                 "subnet_name",
                                 "degree",
                                 "pval",
                                 "fdr",
                                 "annotation"
                               ),
                               draw_heatmap = TRUE,
                               gene_expression_table = NULL,
                               #"~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_log.txt",
                               sample_group_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_sample_description.csv",
                               displayed_groups = c("GSC", "NSC", "NHA"),
                               n_top_genes = 50,
                               importance_field = "score",
                               is_log_gene_expression_table = FALSE,
                               n_top_pathway_to_draw_heatmap = 5,
                               baseline_samples = NULL,
                               annotation_sources = NULL,
                               use_cluster = FALSE,
                               baseline_table = NULL)
{
  #' Extract subnet data
  #' Get comprehensive information for each subnet:
  #'  heatmap,
  #'  expresion values of all genes and significant genes,
  #'   external gene expression values,
  #'    pathway annotations for this subnet. For each pathway in top 5 pathways in this subnet, make heatmap and network
  #'
  #'    @param final_layout data frame of layout from save_layout
  #'    @param subnet_attr
  #'    @param subnet_annotation_by_coverage
  #'    @param combined_subnet_qvalues,
  #'    @param combined_gene_geneset_table,
  #'    @param gene_description_file character name of file describing the biological function of each genes, pulled from RefSeq annotation
  #'    @param qvalue_cutoff numeric the cut off used to define the significant of pathway
  #'    @param cut_off_field = "logFC" ,
  #'    @param sort_fields = c(cut_off_field , "score"),
  #'    @param cut_off_value = 1,
  #'    @param selected_fields character vector of fields to be included in short table
  #' @param outdir  output directory for heatmaps
  #' @param baseline_samples character vector a vector of names of baseline samples for heatmap normalization. The expresion value of each gene for each sample will be divided to the averge value of baseline_samples.
  #'If NULL, then average of all samples will be used as baseline to normalize gene expression values
  
  #' @param gene_expression_table character file name of table containing gene expression in format: SYMBOL, sample1, sample2, .... Rows are gene expression values for each gene
  #' @param sample_group_table character Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ...
  #' @param displayed_groups a character vector: groups of samples to be displayed in heatmap of gene expression among different groups
  #' @param n_top_genes integer Number of top genes to display in heatmaps
  #' @param importance_field character Field to select top genes
  
  #' @param use_cluster bool use cluster for parallezing this fucntion
  #' @return None Just output varyious files  in outdir folder
  dir.create(outdir)
  #outdir = paste(outdir, "/", subnet, sep = "")
  #dir.create(outdir)
  tryCatch({
    if (is.null(final_layout$is_external[1]))
    {
      final_layout = get_external_compartment_field(final_layout)
    }
    if (is.null(final_layout$description[1]))
    {
      final_layout = get_gene_description (final_layout, subnets_attr, gene_description_file)
      
    }
    #table
    n_genes = subnets_attr[subnets_attr$graph_name == subnet, "n_genes"]
    
    if (subnet != "1" )
      # not root subnet
    {
      tryCatch({
        # subnet annotation
        subnet_annotation_table = subnet_annotation_by_coverage[, c("annotation_id", "annotation", subnet)]
        subnet_annotation_table = subnet_annotation_table[subnet_annotation_table [, subnet] >
                                                            0,]
        if (nrow(subnet_annotation_table) > 0)
        {
        subnet_annotation_qvalues = combined_subnet_qvalues[combined_subnet_qvalues$subnet_name ==
                                                              subnet, ]
        subnet_annotation_qvalues = data.frame(t(subnet_annotation_qvalues[, 4:ncol(subnet_annotation_qvalues)]))
        colnames(subnet_annotation_qvalues)[1] = "annotation_q_value"
        subnet_annotation_qvalues$annotation_id = rownames(subnet_annotation_qvalues)
        subnet_annotation_table = fast_merge(
          subnet_annotation_table,
          subnet_annotation_qvalues,
          by.x = "annotation_id",
          by.y = "annotation_id",
          all.x = TRUE
        )
        colnames(subnet_annotation_table)[3] = "pathway_coverage"
        subnet_annotation_table = subnet_annotation_table[order(-subnet_annotation_table$pathway_coverage),]
        subnet_annotation_table$n_genes_in_pathway =  subnet_annotation_table$pathway_coverage * n_genes
        }
        else{
          subnet_annotation_table = NULL
        }
      },
      error = function(err)
      {
        print(paste("3663 network_layout_Son", err))
      })
    }
    genes = final_layout[substr(final_layout$subnet_name, 1, nchar(subnet)) ==
                           subnet,]
    gene_names = as.character(genes$node_name)
    if (subnet != "1" & ! is.null(subnet_annotation_table))
      # not root subnet
    {
      pathway_genes = c()
      for (i in c(1:nrow(subnet_annotation_table)))
      {
        annotation_id = subnet_annotation_table[i, "annotation_id"]
        annotation_genes = combined_gene_geneset_table[, c("SYMBOL", annotation_id)]
        annotation_genes[is.na(annotation_genes)] = 0
        annotation_genes = annotation_genes[annotation_genes[, annotation_id] ==
                                              1, "SYMBOL"]  # all genes in a GO/MSigDB pathway
        share_genes = intersect(gene_names, annotation_genes)  # genes that in both subnet and pathway
        share_genes_string = paste(share_genes, collapse = ", ")
        pathway_genes = c(pathway_genes, share_genes_string)
        if ((i <= n_top_pathway_to_draw_heatmap) &
            draw_heatmap == TRUE)
        {
          annotation = subnet_annotation_table[i, "annotation"]
          
          title = paste("Pathway_", annotation_id, annotation)
          title = gsub("\\(", "_", title)
          title = gsub("\\)", "_", title)
          title = gsub("\\,", "_", title)
          title = gsub("\\/", "_", title)
          pathways_dir = paste0(outdir,"/pathways")
          dir.create(pathways_dir)
          pathway_dir = paste0(pathways_dir, "/", gsub(":", "", annotation_id))
          dir.create(pathway_dir)
          extract_pathway_data(
            final_layout = final_layout,
            annotation_sources = annotation_sources,
            pathway_Id = annotation_id,
            outdir = paste0(pathway_dir, "/global"),
            gene_expression_table = gene_expression_table,
            is_log_gene_expression_table = is_log_gene_expression_table,
            sample_group_table = sample_group_table,
            displayed_groups = displayed_groups,
            n_top_genes = n_top_genes,
            importance_field = importance_field,
            baseline_samples = baseline_samples,
            use_cluster = use_cluster,
            main_title = title,
            baseline_table = baseline_table
          )
          
          draw_heatmap_and_network_for_gene_list(
            gene_list = share_genes,
            outdir = paste0(pathway_dir, "/in_subnet"),
            layout = final_layout,
            gene_expression_table = gene_expression_table,
            is_log_gene_expression_table = is_log_gene_expression_table,
            displayed_groups = displayed_groups,
            n_top_genes = n_top_genes,
            importance_field = importance_field,
            main_title = paste(title,
                               "in subnet",
                               subnet),
            sample_group_table = sample_group_table,
            baseline_samples = baseline_samples,
            use_cluster = use_cluster,
            baseline_table = baseline_table
          )
        }
      }
      subnet_annotation_table$genes = pathway_genes
      filename = paste(outdir, "/subnet_annotation_table.txt", sep = "")
      write.table(
        subnet_annotation_table,
        file = filename,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
    }else
    {
      filename = paste(outdir, "/subnet_annotation_table.txt", sep = "")
      sys(paste("touch", filename))
    }
    
    # order by subnet names, then by sort _field
    
    genes_external = genes[genes$is_external == TRUE, ]
    sign_genes = genes[genes$subnet_qvalue < qvalue_cutoff &
                         abs(genes[, cut_off_field]) > cut_off_value,]
    sign_genes_external = sign_genes[sign_genes$is_external == TRUE, ]
    
    for (field in sort_fields)
    {
      genes_external = genes_external[order(-genes_external[, field]), ]
      sign_genes = sign_genes[order(-sign_genes[, field]), ]
      sign_genes_external = sign_genes_external[order(-sign_genes_external[, field]), ]
      
      # write table
      
      genes = genes[order(-genes[, field]), ]
      filename = paste(outdir, "/genes_sorted_by_", field, ".txt", sep = "")
      write.table(
        genes,
        file = filename,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
      
      print(2576)
      print(cut_off_field)
      print(selected_fields)
      
      print(colnames(genes))
      genes_short = genes[, selected_fields]
      filename = paste(outdir, "/genes_sorted_by_", field, "_short.txt", sep =
                         "")
      write.table(
        genes_short,
        file = filename,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
      
      
      genes = genes[order(genes[, "subnet_name"],-genes[, field]), ]
      filename = paste(outdir,
                       "/genes_sorted_by_subnet_name_and_",
                       field,
                       ".txt",
                       sep = "")
      write.table(
        genes,
        file = filename,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
      
      genes_short = genes[, selected_fields]
      filename = paste(outdir,
                       "/genes_sorted_by_subnet_name_and_",
                       field,
                       "_short.txt",
                       sep = "")
      write.table(
        genes_short,
        file = filename,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
      
      
      
      
      genes_external = genes_external[, selected_fields]
      filename = paste(outdir,
                       "/genes_external_sorted_by_",
                       field,
                       ".txt",
                       sep = "")
      write.table(
        genes_external,
        file = filename,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
      
      sign_genes = sign_genes[, selected_fields]
      filename = paste(outdir,
                       "/significant_genes_sorted_by_",
                       field,
                       ".txt",
                       sep = "")
      write.table(
        sign_genes,
        file = filename,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
      
      sign_genes_external = sign_genes_external[, selected_fields]
      filename = paste(outdir,
                       "/significant_genes_external_sorted_by_",
                       field,
                       ".txt",
                       sep = "")
      write.table(
        sign_genes_external,
        file = filename,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
      
      
    }
    
    
    #normal
    
    # extracellular
    
    # heatmap
    print("4129 network layout Son")
    print(paste("draw heatmap:", draw_heatmap))
    if (draw_heatmap)
    {
      tryCatch({
        get_subnet_heatmap(
          subnet_name = subnet,
          outdir = paste0(outdir,"/subnet_heatmap"),
          layout = final_layout,
          gene_expression_table = gene_expression_table,
          sample_group_table = sample_group_table,
          displayed_groups = displayed_groups,
          n_top_genes = n_top_genes,
          importance_field = importance_field,
          is_log_gene_expression_table = is_log_gene_expression_table,
          baseline_samples = baseline_samples,
          use_cluster = use_cluster,
          baseline_table = baseline_table
        )
      })
    }
    
    #network
    
    print(paste("subnet data extraction for",   subnet, " is  done"))
    
    print(" 3843, done with subnet extraction")
    
  },
  error = function(err)
  {
    print(err)
    
    inputs = list(
      final_layout = final_layout,
      subnet = subnet,
      subnets_attr = subnets_attr,
      subnet_annotation_by_coverage = subnet_annotation_by_coverage,
      combined_subnet_qvalues = combined_subnet_qvalues,
      combined_gene_geneset_table = combined_gene_geneset_table,
      gene_description_file = gene_description_file,
      qvalue_cutoff = qvalue_cutoff,
      outdir = outdir,
      cut_off_field = cut_off_field,
      sort_fields = sort_fields,
      cut_off_value = cut_off_value,
      selected_fields = selected_fields,
      draw_heatmap = draw_heatmap,
      gene_expression_table = gene_expression_table,
      sample_group_table = sample_group_table,
      displayed_groups = displayed_groups,
      n_top_genes = n_top_genes,
      importance_field =  importance_field,
      is_log_gene_expression_table = is_log_gene_expression_table,
      n_top_pathway_to_draw_heatmap = n_top_pathway_to_draw_heatmap,
      baseline_samples = baseline_samples,
      annotation_sources = annotation_sources,
      use_cluster = use_cluster,
      baseline_table = baseline_table
      
    )
    rdsfile = paste0(outdir, "/extract_subnet_data.RDS")
    image_file = paste0(outdir, "/extract_subnet_data.Rdata")
    saveRDS(inputs, file = rdsfile)
    save.image(file = image_file)
    print(paste(
      "error happened at subnet:",
      subnet,
      " data is saved at",
      rdsfile
    ))
  })
  
}




extract_pathway_data = function(final_layout,
                                annotation_sources,
                                pathway_Id = "M5930",
                                outdir = "../data",
                                gene_expression_table = "../data/Math_dataset/Net_Zene/cpm_table.csv",
                                is_log_gene_expression_table = FALSE,
                                sample_group_table = "../data/Math_dataset/sample_groups.csv",
                                displayed_groups =  NULL,
                                n_top_genes = 50,
                                importance_field = "score",
                                baseline_samples = NULL,
                                use_cluster = FALSE,
                                main_title = main_title,
                                baseline_table = NULL)
#' @param baseline_table character name of table describing for each sample which other samples should be used as baseline to rescale for heatmap drawing.
{
  #' pathway is predefined term from GO, MsigDB ID. Extract data for all genes that belong to the pathways
  genes = get_genes_in_pathway(pathway_Id, annotation_sources)
  pathway_data = final_layout[final_layout$node_name %in% genes, ]
  outdir = gsub(":", "", outdir)
  dir.create(outdir)
  outfile = paste(outdir, "/pw_data.txt", sep = "")
  write.table(
    pathway_data,
    file = outfile,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  
  draw_heatmap_and_network_for_gene_list(
    gene_list = genes,
    outdir = outdir,
    layout = final_layout,
    gene_expression_table = gene_expression_table,
    is_log_gene_expression_table = is_log_gene_expression_table,
    sample_group_table = sample_group_table,
    displayed_groups =  displayed_groups ,
    n_top_genes = n_top_genes,
    importance_field = importance_field,
    main_title = main_title,
    baseline_samples = baseline_samples,
    use_cluster = use_cluster,
    baseline_table = baseline_table
  )
  return(pathway_data)
}


# get all gene name in a pathway ( a term in GO or MSigDB) ID
get_genes_in_pathway = function(pathway_ID = "M5930", annotation_sources)
{
  # go = list(gene_geneset_table =go_gene_geneset_table, geneset_description= go_geneset_description, definition_field= "TERM", annotation_id_field ="GOID", annotation_group_field="ONTOLOGY")
  # MSigDB = list(gene_geneset_table =MSigDB_df, geneset_description= MSigDB_description, definition_field= "DESCRIPTION_BRIEF", annotation_id_field ="SYSTEMATIC_NAME", annotation_group_field=NULL)
  # annotation_sources = list(go, MSigDB)
  genes = NULL
  for (i in c(1:length(annotation_sources)))
  {
    #print(i)
    source = annotation_sources[[i]]
    gene_geneset_table = source$gene_geneset_table
    all_pathways = colnames(gene_geneset_table)
    if (pathway_ID %in% all_pathways)
    {
      genes = gene_geneset_table[gene_geneset_table[, pathway_ID] == 1, "SYMBOL"]
      break()
    }
    
  }
  return(genes)
}



full_pipeline = function(network_file = "../data/nha_gsc_net_nbhscore.graphml" ,
                         node_attributes = NULL,
                         graph_format = "graphml",
                         clustering_method = "louvain",
                         max_cluster_size = 100 ,
                         ndim = 2,
                         layout_method = "fr",
                         space_ratios = c(1, 6),
                         node_score_field = "diff",
                         enrichment_score_field = c("node_smooth_score"),
                         sampling_size = 1000,
                         beta = 0.5,
                         n_top_annotation = 3,
                         importance_field = "score",
                         n_top_genes = 3,
                         outdir = "../data",
                         selected_subnets = NULL,
                         #c("1.7", "1.5.5.3"),
                         selected_genes = NULL,
                         #c("NKX6-2","ASCL1","MYCN","BASP1"),
                         gene_description_file = paste0(data_path, "/all_gene_summary.csv"),
                         return_result = FALSE,
                         draw_heatmap = TRUE,
                         gene_expression_table = NULL,
                         # "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_log.txt",
                         sample_group_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_sample_description.csv",
                         displayed_groups = c("GSC", "NSC", "NHA"),
                         n_top_genes_heatmap = 50,
                         selected_fields = c(
                           "node_name",
                           "score",
                           node_score_field,
                           "description",
                           "summary",
                           "subnet_name",
                           "degree",
                           "pvalue",
                           "fdr",
                           "annotation"
                         ),
                         is_log_gene_expression_table = FALSE,
                         do_extract_subnet_data = FALSE,
                         baseline_samples = NULL,
                         use_cluster = FALSE,
                         export_top_master_genes_data = FALSE,
                         baseline_table = NULL,
                         network_dir="network")


{
  dir.create(outdir)
  job_log_dir <<-
    paste0(outdir, "/jobs")  # Folder containing job files to submit to computer cluster. Global variable
  dir.create(job_log_dir)
  print(paste("4233 network_layout_Son job_log_dir", job_log_dir))
  #outdir = paste(outdir, "/layout", sep = "")
  outdir = paste(outdir, "/",network_dir, sep = "")
  
  dir.create(outdir)
  print(node_attributes[1:5,])
  node_attributes = node_attributes[!duplicated(node_attributes[, 1]), ]  #remove duplicated node_attributes rows
  #print("4274 node attributes:")
  #print(node_attributes[1:5, ])
  layout_dir = paste0(outdir, "/layout")
  dir.create(layout_dir)
  for (i in c(1:length(space_ratios)))
  {
    space_ratio = space_ratios[i]
    print(paste("space ratio:", space_ratio))

    layout_dir_space_ratio = paste(layout_dir, "/", space_ratio, sep = "")
    pipeline_out = get_network_layout_pipeline(
      network_file = network_file,
      node_attributes = node_attributes,
      graph_format = graph_format,
      clustering_method = clustering_method,
      max_cluster_size = max_cluster_size,
      ndim = ndim,
      layout_method = layout_method,
      space_ratio = space_ratio,
      node_score_field = node_score_field,
      enrichment_score_field = enrichment_score_field,
      sampling_size = sampling_size,
      beta = beta,
      n_top_annotation = n_top_annotation,
      importance_field = importance_field,
      n_top_genes = n_top_genes,
      outdir = outdir
    )
    if (i != 1)
    {
      draw_heatmap = draw_heatmap
      do_extract_subnet_data = do_extract_subnet_data
      save_network_layout_only = FALSE
    }
    else
    {
      draw_heatmap = FALSE
      do_extract_subnet_data = FALSE
      save_network_layout_only = TRUE
    }
    
    
    out = save_layout(
      outdir = outdir,
      gene_description_file = gene_description_file,
      pipeline_out = pipeline_out,
      selected_subnets = selected_subnets,
      # immune subnet,
      selected_genes = selected_genes,
      qvalue_cutoff = 0.05,
      sort_fields = c(node_score_field, "score"),
      cut_off_field = node_score_field,
      cut_off_value = 1,
      selected_fields = selected_fields,
      draw_heatmap = draw_heatmap,
      gene_expression_table = gene_expression_table,
      sample_group_table = sample_group_table,
      displayed_groups = displayed_groups,
      n_top_genes = n_top_genes_heatmap,
      importance_field = importance_field,
      is_log_gene_expression_table = is_log_gene_expression_table,
      do_extract_subnet_data = do_extract_subnet_data,
      baseline_samples = baseline_samples,
      use_cluster = use_cluster,
      save_network_layout_only = save_network_layout_only,
      export_top_master_genes_data = export_top_master_genes_data,
      baseline_table = baseline_table,
      layout_dir = layout_dir_space_ratio
    )
  }
  
  
}




full_pipeline_one_space_ratio = function(network_file = "../data/nha_gsc_net_nbhscore.graphml" ,
                                         node_attributes = NULL,
                                         graph_format = "graphml",
                                         clustering_method = "louvain",
                                         max_cluster_size = 100 ,
                                         ndim = 2,
                                         layout_method = "fr",
                                         space_ratios = 1,
                                         node_score_field = "diff",
                                         enrichment_score_field = c("node_smooth_score"),
                                         sampling_size = 1000,
                                         beta = 0.5,
                                         n_top_annotation = 3,
                                         importance_field = "score",
                                         n_top_genes = 3,
                                         outdir = "../data",
                                         selected_subnets = NULL,
                                         #c("1.7", "1.5.5.3"),
                                         selected_genes = NULL,
                                         #c("NKX6-2","ASCL1","MYCN","BASP1"),
                                         gene_description_file = paste0(data_path, "/all_gene_summary.csv"),
                                         return_result = FALSE,
                                         draw_heatmap = TRUE,
                                         gene_expression_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_log.txt",
                                         sample_group_table = "~/Desktop/SonData/BlenderExercises/vtk_test/data/GSC_NSC_NHA_normalized_gene_expression_sample_description.csv",
                                         displayed_groups = c("GSC", "NSC", "NHA"),
                                         n_top_genes_heatmap = 50,
                                         selected_fields = c(
                                           "node_name",
                                           "score",
                                           node_score_field,
                                           "description",
                                           "summary",
                                           "subnet_name",
                                           "degree",
                                           "pvalue",
                                           "fdr",
                                           "annotation"
                                         ),
                                         is_log_gene_expression_table = FALSE,
                                         do_extract_subnet_data = FALSE)


{
  outdir_space_ratio = paste(outdir, "/", space_ratio, sep = "")
  pipeline_out = get_network_layout_pipeline(
    network_file = network_file,
    node_attributes = node_attributes,
    graph_format = graph_format,
    clustering_method = clustering_method,
    max_cluster_size = max_cluster_size,
    ndim = ndim,
    layout_method = layout_method,
    space_ratio = space_ratio,
    node_score_field = node_score_field,
    enrichment_score_field = enrichment_score_field,
    sampling_size = sampling_size,
    beta = beta,
    n_top_annotation = n_top_annotation,
    importance_field = importance_field,
    n_top_genes = n_top_genes,
    outdir = outdir_space_ratio
  )
  
  out = save_layout(
    outdir = outdir_space_ratio,
    gene_description_file = gene_description_file,
    pipeline_out = pipeline_out,
    selected_subnets = selected_subnets,
    # immune subnet,
    selected_genes = selected_genes,
    qvalue_cutoff = 0.05,
    sort_fields = c(node_score_field, "score"),
    cut_off_field = node_score_field,
    cut_off_value = 1,
    selected_fields = selected_fields,
    draw_heatmap = draw_heatmap,
    gene_expression_table = gene_expression_table,
    sample_group_table = sample_group_table,
    displayed_groups = displayed_groups,
    n_top_genes = n_top_genes_heatmap,
    importance_field = importance_field,
    is_log_gene_expression_table = is_log_gene_expression_table,
    do_extract_subnet_data = do_extract_subnet_data
  )
  
  
}

#Gene expression table in format ( SYMBOL, SAMPLE1, SAMPLE2, etc, row are expression value for each gene)
get_subnet_heatmap = function(subnet_name = "1.2.1",
                              outdir = "../data/Math_dataset/Net_Zene/Fibroblast2ESC/layout/1",
                              layout = "../data/Math_dataset/Net_Zene/Fibroblast2ESC/layout/1/layout.txt",
                              gene_expression_table = "../data/Math_dataset/Net_Zene/cpm_table.csv",
                              is_log_gene_expression_table = FALSE,
                              sample_group_table = "../data/Math_dataset/sample_groups.csv",
                              displayed_groups =  NULL,
                              #c("Fibroblast", "ESC", "iPSC", "GL261_Control"),
                              n_top_genes = 50,
                              importance_field = "score",
                              baseline_samples = NULL,
                              use_cluster = FALSE,
                              baseline_table = NULL)
{
  #' Get heatmaps for a subnet
  
  #' Drawing 3 type of heatmap : all genes, significant genes, top genes both as static image or html interactive image
  #'Gene expression table in format ( SYMBOL, SAMPLE1, SAMPLE2, etc, row are expression value for each gene)
  
  #' @param subnet_name character Name of subnet to display
  #' @param outdir  output directory for heatmaps
  #' @param layout character name of network layout file
  #' @param gene_expression_table character file name of table containing gene expression in format: SYMBOL, sample1, sample2, .... Rows are gene expression values for each gene
  #' @param is_log_gene_expression_table Boolean TRUE/FALSE if values in gene expression tables are in log format. If not, then function will automatically transform into log format
  #' @param sample_group_table character Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ...
  #' @param displayed_groups a character vector: groups of samples to be displayed in heatmap of gene expression among different groups
  #' @param n_top_genes integer Number of top genes to display in heatmaps
  #' @param importance_field character Field to select top genes
  #' @param baseline_samples character vector a vector of names of baseline samples for heatmap normalization. The expresion value of each gene for each sample will be divided to the averge value of baseline_samples.
  #'If NULL, then average of all samples will be used as baseline to normalize gene expression values
  #' @param use_cluster bool use cluster for parallezing this fucntion
  #' @param baseline_table character name of table describing for each sample which other samples should be used as baseline to rescale for heatmap drawing.
  #' @return None Just draw heatmaps in outdir folder in html and png format
  
  get_heatmaps(
    outdir = outdir,
    layout = layout,
    gene_expression_table = gene_expression_table,
    is_log_gene_expression_table = is_log_gene_expression_table,
    sample_group_table = sample_group_table,
    displayed_groups =  displayed_groups ,
    n_top_genes = n_top_genes,
    importance_field = importance_field ,
    main_title = paste("Heatmap of subnetwork", subnet_name),
    gene_selection_method = "subnet_name",
    subnet_name = subnet_name,
    gene_list = NULL,
    baseline_samples = baseline_samples,
    use_cluster = use_cluster,
    baseline_table = baseline_table
  )
  
}

draw_heatmap_and_network_for_gene_list = function(gene_list = c("POU5F1", "MYCN", "KLF4", "SOX2"),
                                                  outdir = "../data/Math_dataset/Net_Zene/Fibroblast2ESC/layout/1",
                                                  layout = "../data/Math_dataset/Net_Zene/Fibroblast2ESC/layout/1/layout.txt",
                                                  gene_expression_table = "../data/Math_dataset/Net_Zene/cpm_table.csv",
                                                  is_log_gene_expression_table = FALSE,
                                                  sample_group_table = "../data/Math_dataset/sample_groups.csv",
                                                  displayed_groups =  NULL,
                                                  n_top_genes = 50,
                                                  importance_field = "score",
                                                  main_title = "Heatmap of subnetwork",
                                                  baseline_samples = NULL,
                                                  use_cluster = FALSE,
                                                  baseline_table = NULL)
#c("Fibroblast", "ESC", "iPSC", "GL261_Control")
#' Draw heatmap and network for a gene list
#' Given a gene list, draw heatmap of these genes and the network containing only these genes.
#' @param gene_list character vector of genes to display
#' @param outdir  output directory for result
#' @param layout character name of network layout file
#' @param gene_expression_table character file name of table containing gene expression in format: SYMBOL, sample1, sample2, .... Rows are gene expression values for each gene
#' @param is_log_gene_expression_table Boolean TRUE/FALSE if values in gene expression tables are in log format. If not, then function will automatically transform into log format
#' @param sample_group_table character Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ...
#' @param displayed_groups a character vector: groups of samples to be displayed in heatmap of gene expression among different groups
#' @param n_top_genes integer Number of top genes to display in heatmaps
#' @param importance_field character Field to select top genes
#' @param main_title character prefix of title of heatmap
#' @param baseline_samples character vector a vector of names of baseline samples for heatmap normalization. The expresion value of each gene for each sample will be divided to the averge value of baseline_samples.
#'If NULL, then average of all samples will be used as baseline to normalize gene expression values
#' @param job_log_dir character folder for heatmap job scripts
#' @param use_cluster bool use cluster for parallezing this fucntion
#' @param baseline_table character name of table describing for each sample which other samples should be used as baseline to rescale for heatmap drawing.
#' @return None Just draw heatmaps in outdir folder in html and png format

{
  print(gene_list)
  get_heatmaps(
    outdir = outdir,
    layout = layout ,
    gene_expression_table = gene_expression_table,
    is_log_gene_expression_table = is_log_gene_expression_table,
    sample_group_table = sample_group_table,
    displayed_groups =  displayed_groups,
    n_top_genes = n_top_genes,
    importance_field = importance_field,
    main_title = main_title ,
    gene_selection_method = "gene_list",
    gene_list = gene_list,
    baseline_samples = baseline_samples,
    use_cluster = use_cluster,
    baseline_table = baseline_table
  )
  
}


get_heatmaps = function(outdir = "../data/Math_dataset/Net_Zene/Fibroblast2ESC/network_layout_result/1",
                        layout = "../data/Math_dataset/Net_Zene/Fibroblast2ESC/network_layout_result/1/layout.txt",
                        gene_expression_table = "../data/Math_dataset/Net_Zene/cpm_table.csv",
                        is_log_gene_expression_table = FALSE,
                        sample_group_table = "../data/Math_dataset/sample_groups.csv",
                        displayed_groups =  NULL,
                        n_top_genes = 50,
                        importance_field = "score",
                        main_title = "Heatmap of subnetwork",
                        gene_selection_method = "subnet_name",
                        # c("subnet_name", "gene_list) How genes are extracted to draw heatmap. It can be by subnet names (option "subnet_name", or by providing gene list ("gene_list")
                        subnet_name = NULL,
                        gene_list = NULL,
                        baseline_samples = NULL,
                        use_cluster = FALSE,
                        baseline_table = NULL, # "../data/Math_dataset/baselines.csv"
                        show_subnet_name_in_heatmap=TRUE)

{
  #' Get heatmaps
  
  #' Drawing 3 type of heatmap : all genes, significant genes, top genes both as static image or html interactive image for gene list or a subnet depending on the gene_selection_method
  #'Gene expression table in format ( SYMBOL, SAMPLE1, SAMPLE2, etc, row are expression value for each gene)
  
  #' @param baseline_samples character vector a vector of names of baseline samples for heatmap normalization. The expresion value of each gene for each sample will be divided to the averge value of baseline_samples.
  #'If NULL, then average of all samples will be used as baseline to normalize gene expression values. It can also be a string of comma separated list of samples like "Sample1, Sample2, Sample3".
  #' @param outdir  output directory for heatmaps
  #' @param layout character name of network layout file
  #' @param gene_expression_table character file name of table containing gene expression in format: SYMBOL, sample1, sample2, .... Rows are gene expression values for each gene
  #' @param is_log_gene_expression_table Boolean TRUE/FALSE if values in gene expression tables are in log format. If not, then function will automatically transform into log format
  #' @param sample_group_table character Table describing samples and factors(groups). Contains  columns: Sample, Factor1, Factor2, ...
  #' @param displayed_groups a character vector: groups of samples to be displayed in heatmap of gene expression among different groups
  #' @param n_top_genes integer Number of top genes to display in heatmaps
  #' @param importance_field character Field to select top genes
  #' @param main_title character prefix of title of heatmap  such as"Heatmap of subnetwork",
  #' @param gene_selection_method character  c("subnet_name", "gene_list) How genes are extracted to draw heatmap. It can be by subnet names (option "subnet_name", or by providing gene list ("gene_list")
  #' @param subnet_name character Name of subnet to display
  #' @param gene_list character vector of genes to display
  #' @param baseline_samples character vector of samples to make baseline to normalize gene expression values
  #' @param use_cluster bool use cluster for parallezing this fucntion
  #' @param baseline_table character name of table describing for each sample which other samples should be used as baseline to rescale for heatmap drawing.
  #' Contain two colulmn: "SampleID, Baseline_samples", where in the column Baseline_samples is the comma separated list of basseline samples for corresponding SampleID. This table takes priority over baseline_samples parameter
  #' @param show_subnet_name_in_heatmap bool show subnet name after gene name in heatmap
  #' @return None Just draw heatmaps in outdir folder in html and png format
  #' 
 print(paste("4619 network_layout_Son baseline_table:", baseline_table))
   inputs = list(
    subnet_name = subnet_name ,
    outdir = outdir,
    layout = layout,
    gene_expression_table = gene_expression_table,
    is_log_gene_expression_table = is_log_gene_expression_table,
    sample_group_table = sample_group_table,
    displayed_groups =  displayed_groups,
    n_top_genes = n_top_genes,
    importance_field = importance_field,
    main_title = main_title,
    gene_selection_method = gene_selection_method,
    gene_list = gene_list,
    baseline_samples = baseline_samples,
    use_cluster = use_cluster
  )
  if (class(baseline_samples) == "character")
  {
    baseline_samples = unlist(strsplit(baseline_samples, ","))
    baseline_samples = trimws(baseline_samples)
  }
  main_title = gsub(":", "_", main_title) # Change the name of the title , replace : with _ to make it work for Windows if the main title contains GO:xxx.
  #saveRDS(inputs, file=paste0(outdir,"/get_heatmap.RDS"))
  #save.image(file = paste0(outdir,"/get_heatmap.Rdata"))
  print("4066 getsubnetheatmap network_layout_Son_053018.R")
  print(paste("heatmap title:", main_title))
  dir.create(outdir)
  # suboutdir = gsub(" ", "_", main_title)
  # outdir = paste0(outdir, "/", suboutdir)
  # dir.create(outdir)
  # import imput files
  if (class(layout) == "character")
    # layout as a text file:
  {
    layout <- import_file(layout)
  }
  duplicated_genes = duplicated(layout$node_name)
  layout = layout[!duplicated_genes,]  # remove duplicated genes
  
  if (gene_selection_method == "subnet_name")
    # Select layout . Can be either by subnet_name or gene list
  {
    if (!is.null(subnet_name))
    {
      layout_selected = data.frame(layout[substr(layout$subnet_name, 1, nchar(subnet_name)) ==
                                            subnet_name, ])
      
      
    } else
    {
      stop("4456, network_layout_Son subnet name has NOT been provided")
    }
  } else
  {
    # using gene list
    if (!is.null(gene_list))
    {
      gene_list = toupper(gene_list)
      layout$node_name = toupper(layout$node_name)
      layout_selected = data.frame(layout[layout$node_name %in% gene_list, ])
    } else{
      stop("4465 network_layout_Son gene list has NOT been provided")
    }
    
    
  }
  
  
  gene_expression_table <-
    import_file(gene_expression_table)
  colnames(gene_expression_table)[1] = "SYMBOL"
  gene_expression_table$SYMBOL = toupper(gene_expression_table$SYMBOL)
  sample_list = list()
  
  sample_order_original = NULL
  if (class(sample_group_table) == "character")
    # layout as a text file:
  {
    sample_group_table <- import_file(sample_group_table)
    sample_order_original =  unlist(sample_group_table[ ,1])
    if (ncol(sample_group_table) == 2)
    {
      colnames(sample_group_table) = c("Sample", "Group")
    } else
    {
      groups = apply(sample_group_table[, 2:ncol(sample_group_table)], 1, function(row)
        paste(row, collapse = "_"))
      sample_group_table = data.frame(Sample = sample_group_table[, 1], Group = groups)
      colnames(sample_group_table) = c("Sample", "Group")
    }
    
  }
  selected_samples = intersect(sample_group_table$Sample,
                               colnames(gene_expression_table))
  sample_group_table = sample_group_table[sample_group_table$Sample %in% selected_samples, ]
  gene_expression_table = gene_expression_table[, c("SYMBOL", selected_samples)]
  
  if (is.null(displayed_groups))
  {
    displayed_groups = unique(sample_group_table$Group)
  }
  
  
  #if (!is.null(displayed_groups))
  #{
  for (group in displayed_groups)
  {
    samples = unlist(sample_group_table[sample_group_table[, 2] == group, 1])
    sample_list[[group]] = samples
  }
  
  sample_list_df = data.frame(Sample = as.character(unlist(sample_list)))
  if (nrow(sample_list_df) == 0)
  {
    print(sample_list_df)
    stop("4158 network_layout_son, sample_list_df has zero rows")
  }
  
  if (nrow(sample_group_table) == 0)
  {
    print(sample_group_table)
    stop("4158 network_layout_son, sample_group_table has zero rows")
  }
  
  sample_list_df = fast_merge(sample_list_df,
                              sample_group_table,
                              by.x = "Sample",
                              by.y = "Sample")
  sample_list_new = unlist(sample_list)
  gene_expression_displayed_groups = gene_expression_table[c("SYMBOL", sample_list_new)]
  
  
  layout_selected = fast_merge(
    layout_selected,
    gene_expression_displayed_groups,
    by.x = "node_name",
    by.y = "SYMBOL",
    all.x = FALSE,
    all.y = FALSE
  )
  
  if (!is.null(sample_order_original))
  {# now need to order sample_list_new accoring to sample_order_original to keep order when drawing heatmap. 
  sample_list_new = sample_order_original[sample_order_original %in% sample_list_new]
  }
  
  heatmap_data = layout_selected[, sample_list_new]
  
  
  if (show_subnet_name_in_heatmap==TRUE)
  {
    rownames(heatmap_data) = paste(layout_selected$node_name,layout_selected$subnet_name) 
  }else
  {
    rownames(heatmap_data) = layout_selected$node_name
  }
  
  heatmap_data[is.na(heatmap_data)] = 0
  
  heatmap_rowsum = rowSums(heatmap_data)
  heatmap_zerorow = heatmap_rowsum[heatmap_rowsum == 0]
  zero_genes = names(heatmap_zerorow)  # genes with gene expression equal zero across all samples. Exclude for downstream analysis.
  gene_names = rownames(heatmap_data)
  filtered_genes = gene_names[!gene_names %in% zero_genes]
  print("4760 network_layout_Son_053018.R filtered gene names:")
  print(filtered_genes)
  heatmap_data = heatmap_data[filtered_genes,]
  
  layout_selected = layout_selected[paste(layout_selected$node_name, layout_selected$subnet_name) %in% filtered_genes, ] # only include filtered genes
  ngene = nrow(layout_selected)
 

  
  
  
  
  
  
  # Normalizing data
  if (!is_log_gene_expression_table)
  {
    #Finding smallest value to add before taking log
    expression_unlist = unlist(heatmap_data)
    expression_unlist = as.numeric(expression_unlist)
    expression_unlist = sort(expression_unlist)
    for (i in c(1:length(expression_unlist)))
    {
      min_expression_value = expression_unlist[i]
      if (min_expression_value > 0) # smallest postive number
      {
        break
      }
    }
    # Add this smallest value to heatmap_data so that  zero value would not screw up the average result.
    small_value = 0.5* min_expression_value 
    heatmap_data = heatmap_data + small_value # so that avoid error when taking log of 0 (cpm)
    heatmap_data = log2(heatmap_data)
  }
  heatmap_data = as.data.frame(heatmap_data)
  
  #########################################baseline normalization
  if (!is.null(baseline_table))
  {
    baseline_table_df = fread(baseline_table, data.table = FALSE)
    colnames(baseline_table_df)[1:2] = c("Sample", "Baseline_Samples")
    rownames(baseline_table_df) = baseline_table_df$Sample
    heatmap_data_rescale = heatmap_data
    for (sample in colnames(heatmap_data))
    {
      if (sample %in% baseline_table_df$Sample)
      {
        baselines = baseline_table_df[sample, "Baseline_Samples"]
      } else
      {
        print(paste(
          "sample ",
          sample,
          "is not in baseline_table ",
          baseline_table
        ))
        stop()
      }
      baseline_samples = trimws(unlist(strsplit(baselines, ",")))
      baseline_samples = intersect(baseline_samples, colnames(heatmap_data)) # baseline samples should be in both baseline table and heatmap data
      tryCatch({
        heatmap_data_baseline = heatmap_data[, baseline_samples]
      },
      error = function(cond) {
        message(
          "the given baseline samples list contains samples that are not in the gene expression table. check baseline samples and gene expression table"
        )
        message("Here's the original error message:")
        message(cond)
        stop()
      })
      baseline = apply(heatmap_data_baseline, 1, function(x) {
        median(x, na.rm = TRUE)
      })
      heatmap_data_rescale[, sample] = heatmap_data[, sample] - as.numeric(unlist(baseline))
      
      
    }
    heatmap_data = heatmap_data_rescale
  } else
  {
    if (is.null(baseline_samples))
    {
      baseline_samples = colnames(heatmap_data)
    }
    # print(paste(
    #   "4487 network layout Son, baseline_samples",
    #   baseline_samples
    # ))
    # 
    tryCatch({
      heatmap_data_baseline = heatmap_data[, baseline_samples]
    },
    error = function(cond) {
      message(
        "the given baseline samples list contains samples that are not in the gene expression table. check baseline samples and gene expression table"
      )
      message("Here's the original error message:")
      message(cond)
      stop()
    })
    
    #baseline = rowMeans(heatmap_data_baseline, na.rm = TRUE)
    baseline = apply(heatmap_data_baseline, 1, function(x) {
      median(x, na.rm = TRUE)
    })
    #heatmap_data = t(heatmap_data)  # transpose so can normalized by scale function
    #heatmap_data = scale(heatmap_data)
    #heatmap_data = t(heatmap_data)
    if (nrow(heatmap_data) > 1)
    {
      heatmap_data = as.data.frame(sapply(heatmap_data,  function(x) {
        x - baseline
      }))
    } else
    {
      print("4694, only one gene in heatmap data")
      heatmap_data = heatmap_data - as.numeric(unlist(baseline))
      print(heatmap_data)
    }
  }
  #################################################################################
  if (ngene > n_top_genes)
  {
    layout_selected_ordered = layout_selected[order(layout_selected[, importance_field], decreasing = TRUE),]  # order by importance field
    top_genes_layout = layout_selected_ordered[c(1:n_top_genes), c("node_name", "subnet_name")]
    if (show_subnet_name_in_heatmap==TRUE)
    {

      top_genes = paste(top_genes_layout$node_name,  top_genes_layout$subnet_name)
    }else
    {
      top_genes = top_genes_layout$node_name
    }
    heatmap_data_top_genes = heatmap_data[top_genes,]
  }
  
  
  
  
  if (is.element("pvalue", colnames(layout_selected)))
  {
    significant_genes = layout_selected[layout_selected$pvalue < 0.05 & abs(layout_selected$logFC) > 1, "node_name"]
  }else
  {
    significant_genes = layout_selected[ abs(layout_selected$logFC) > 1, "node_name"]
  }
  # if (is.element("fdr", colnames(layout_selected)))
  # {
  #   significant_genes = layout_selected[layout_selected$fdr < 0.05 & abs(layout_selected$logFC) > 0, "node_name"]
  # }
  
  if (show_subnet_name_in_heatmap==TRUE)
  {
    significant_genes_layout = layout_selected[layout_selected$node_name %in% significant_genes, c("node_name", "subnet_name")]
    significant_genes = paste(significant_genes_layout$node_name,  significant_genes_layout$subnet_name)
  }
  
  heatmap_data_significant_genes = heatmap_data[significant_genes,]
  

  
  
  # Data frame with column annotations.
  # if (!is.null(displayed_groups))
  # {
  mat_col <-
    data.frame(group = sample_list_df$Group)
  generated_colors = brewer.pal(length(displayed_groups), "Set1")
  if (length(generated_colors) > length(displayed_groups))
  {
    generated_colors = generated_colors[1:length(displayed_groups)]
  }
  if (length(generated_colors) < length(displayed_groups))
  {
    # pallete does not have enough number of colors for displayed groups. generate pallete using other strategy
    # manually create color range
    myColors = c("green", "yellow", "red")
    
    # expand the color range
    generated_colors = colorRampPalette(myColors)(length(displayed_groups))
 
    
    
  }
  mat_colors = sample(colors(), length(displayed_groups), replace = FALSE)
  print("mat colors:")
  print(mat_colors)
  mat_colors <- list(group= mat_colors) #  list(group = generated_colors)
  
  print("mat colors:")
  print(mat_colors)
  
  names(mat_colors$group) <-
    unique(mat_col$group)
  
  # } else
  # {
  #   mat_col = NULL
  #   mat_colors = NULL
  # }
  #
  
  for (cluster_col in c(FALSE, TRUE))
  {
    cluster_state = ""
    if (cluster_col == TRUE)
    {
      cluster_state = "_clustered_"
    }
    if (n_top_genes < ngene)
      # number of genes in this subset data shoud be less than total number of gene in all gene heatmap before heatmap for data subset been drawed
    {
      display_heatmap(
        data = heatmap_data_top_genes,
        file = paste0("top_genes", cluster_state, ".txt"),
        cluster_row = FALSE,
        cluster_col = cluster_col,
        title = paste(main_title,
                      "top",
                      n_top_genes,
                      "genes"),
        mat_col = mat_col,
        mat_colors = mat_colors,
        outdir = outdir,
        use_cluster = use_cluster
      )
      
    }
    
    display_heatmap(
      data = heatmap_data,
      file = paste0("all_genes", cluster_state, ".txt"),
      cluster_row = TRUE,
      cluster_col = cluster_col,
      title = paste(main_title, " all genes"),
      mat_col = mat_col,
      mat_colors = mat_colors,
      outdir = outdir,
      use_cluster = use_cluster
    )
    
    
    
    display_heatmap(
      data = heatmap_data_significant_genes,
      file = paste0("sign_genes", cluster_state, ".txt"),
      cluster_row = TRUE,
      cluster_col = cluster_col,
      title = paste(main_title, "significantly_changed_genes"),
      mat_col = mat_col,
      mat_colors = mat_colors,
      outdir = outdir,
      use_cluster = use_cluster
    )
    
    
    display_heatmap(
      data = heatmap_data,
      file = paste0("all_genes_subnet_order", cluster_state, ".txt"),
      cluster_row = FALSE,
      cluster_col = cluster_col,
      title = paste(main_title, " all genes with subnet order"),
      mat_col = mat_col,
      mat_colors = mat_colors,
      outdir = outdir,
      use_cluster = use_cluster
    )
    
    
    
    display_heatmap(
      data = heatmap_data_significant_genes,
      file = paste0("sign_genes_subnet_order", cluster_state, ".txt"),
      cluster_row = FALSE,
      cluster_col = cluster_col,
      title = paste(main_title, "significantly_changed_genes subnet order"),
      mat_col = mat_col,
      mat_colors = mat_colors,
      outdir = outdir,
      use_cluster = use_cluster
    )
    
      
    
  }
  
}




option_list <- list(
  make_option(
    "--cmd",
    default = "help",
    help = "run network layout commands . Accepted values:heatmap",
    type = "character"
  ),
  make_option(
    "--heatmap_RDS_file",
    default = NULL,
    help = "file name of heatmap RDS to draw heatmap",
    type = "character"
  )
)

print(paste("4723 network_layout_Son"))


# pipeline_out = get_network_layout_pipeline(ndim = 2, space_ratio = 4, max_cluster_size = 100)
#pipeline_out = readRDS("~/Desktop/SonData/BlenderExercises/vtk_test/data/layout/get_network_layout_pipeline.RDS")
#save_layout(pipeline_out= pipeline_out,selected_subnets =c("1.7", "1.5.5.3")) # immune subnet
# graph = pipeline_out[[1]]
# final_layout = pipeline_out[[2]]
# updated_subnets = pipeline_out[[3]]
# subnets_attr = pipeline_out[[4]]
# subnets_annotations = pipeline_out[[5]]
# sign_subnets =pip
# RUNX1 = check_gene_function(gene_name = "RUNX1", outdir = "../data", layout= final_layout, subnets_attr = subnets_attr)
# Need to do: heatmap for each subnet. Heatmap for each gene annotation/function prediction
#pipeline_out_after_save= full_pipeline(space_ratio=1)
# graph=pipeline_out_after_save$graph
# layout=pipeline_out_after_save$layout
# subnets =pipeline_out_after_save$subnets
# subnets_attr=pipeline_out_after_save$subnets_attr
# annotations=pipeline_out_after_save$annotations
# sign_subnets = pipeline_out_after_save$sign_subnets_attr
# annotation_sources= pipeline_out_after_save$annotation_sources
# combined_annotations = pipeline_out_after_save$combined_annotations
# subnet_top_annotation = pipeline_out_after_save$subnet_top_annotation
# subnet_annotation_by_coverage= pipeline_out_after_save$subnet_annotation_by_coverage
# combined_subnet_qvalues= pipeline_out_after_save$combined_subnet_qvalues
# combined_gene_geneset_table = pipeline_out_after_save$combined_gene_geneset_table
#
# CEBPB=extract_gene_function(gene_name ="CEBPB",
#                       layout= layout,
#                       subnets_attr=subnets_attr,
#                       outdir = paste("../data",sep=""),
#                       subnet_annotation_by_coverage=subnet_annotation_by_coverage,
#                       combined_subnet_qvalues=combined_subnet_qvalues,
#                       combined_gene_geneset_table=combined_gene_geneset_table
# )

#EMT=extract_pathway_data(final_layout=layout, annotation_sources=annotation_sources, pathway_Id="GO:0001837")
# out = readRDS("~/Desktop/SonData/BlenderExercises/vtk_test/data/layout/1/get_network_layout_pipeline_after_saving.RDS")
# graph=out$graph
# layout=out$layout
# subnets =out$subnets
# subnets_attr=out$subnets_attr
# annotations=out$annotations
# sign_subnets = out$sign_subnets_attr
# annotation_sources= out$annotation_sources
# combined_annotations = out$combined_annotations
# subnet_top_annotation = out$subnet_top_annotation
# subnet_annotation_by_coverage= out$subnet_annotation_by_coverage
# combined_subnet_qvalues= out$combined_subnet_qvalues
# combined_gene_geneset_table = out$combined_gene_geneset_table

# STOX2_subnet= extract_subnet_data(final_layout=layout,
#                               subnet="1.4.4",
#                               subnets_attr,
#                               subnet_annotation_by_coverage,
#                               combined_subnet_qvalues,
#                               combined_gene_geneset_table,
#
#                               gene_description_file="../data/all_gene_summary.csv",
#                               qvalue_cutoff=0.05,
#                               outdir="../data/layout",
#                               sort_fields =c("logFC", "score"),
#                               cut_off_field = "logFC",
#                               cut_off_value =1,
#                               selected_fields = c("node_name","score","logFC", "description", "summary","subnet_name","degree", "pval", "pfdr",  "annotation" ),
#                               sample_fields=c("logFC")  # fields containing samples in final layout to draw heatmap,
#
#
# )





# full_pipeline(network_file = "../data/Math_dataset/mouse_embryonic_stem_cell_net.csv",
#               node_attributes= "../data/Math_dataset/Net_Zene/Fibroblast2ESC/nSCORE/Fibroblast2ESC_dea/scores.csv",
#               graph_format = "ncol",
#               clustering_method = "louvain",
#               node_score_field="logFC",
#               outdir="../data/Math_dataset/Net_Zene/Fibroblast2ESC",
#               selected_subnets =NULL,
#               selected_genes=NULL,
#               gene_description_file = "../data/all_gene_summary.csv",
#               return_result=FALSE,
#               draw_heatmap =TRUE,
#               gene_expression_table="../data/Math_dataset/Net_Zene/cpm_table.csv",
#               sample_group_table="../data/Math_dataset/sample_groups.csv",
#               displayed_groups=c("Fibroblast", "ESC", "iPSC", "GL261_Control"),
#               n_top_genes_heatmap=50,
#               selected_fields=c("node_name","score","logFC", "description", "summary","subnet_name","degree", "pvalue", "fdr",  "annotation" ))
#
#
# get_subnet_heatmap(subnet_name="1", outdir="../data/Math_dataset/test_heatmap_row_cluster",
#                             layout="../data/Math_dataset/test_python/layout/1/layout.txt",
#                             gene_expression_table="../data/Math_dataset/Net_Zene/cpm_table.csv",
#                             is_log_gene_expression_table=FALSE,
#                             sample_group_table="../data/Math_dataset/sample_groups.csv",
#                             displayed_groups=c("Fibroblast","ESC","iPSC","GL261_Control"),
#                             n_top_genes=50,
#                             importance_field="score")
