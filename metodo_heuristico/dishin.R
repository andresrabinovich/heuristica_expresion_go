common_ancestors <- function(a, b){
  bla <- intersect(ancestros[[a]], ancestros[[a]])
  return(bla)
}

pot <- function(x, k){
  r <- diag(dim(x)[2])
  for(i in 1:k) r <- r %*% x
  r
}


g <- graph_from_graphnel(p[[1]]$dag)
g <- induced_subgraph(g, sample(V(g), 500))
a <- as.matrix(get.adjacency(g))

caminos <- matrix(0, ncol=vcount(g), nrow=vcount(g))
for (i in 1:100) {
  b <- pot(a, i)
  rownames(b)=colnames(b)
  diag(b) <-rep(0, vcount(g))
  caminos <- caminos + b
}



library("RMySQL")
con<-dbConnect(RMySQL::MySQL(), dbname = "go", username = "root", password = "bambu")
dbSendQuery(con, "SET @t1Id = (SELECT id FROM term WHERE acc='GO:0060255'),  @t2Id = (SELECT id FROM term WHERE acc='GO:0031326');")
dbSendQuery(con, "SET @maxFreq = (SELECT COUNT(*) FROM gene_product);")
dbSendQuery(con, "SET @t1IC = (SELECT -LOG(COUNT(DISTINCT a.gene_product_id)/@maxFreq) as ic FROM graph_path gp INNER JOIN association a ON (gp.term2_id = a.term_id) WHERE gp.term1_id = @t1Id AND a.is_not = 0 AND gp.relationship_type_id IN (SELECT id FROM term WHERE name='part_of' OR name='is_a'));")
dbSendQuery(con, "SET @t2IC = (
  SELECT -LOG(COUNT(DISTINCT a.gene_product_id)/@maxFreq) as ic
  FROM graph_path gp
  INNER JOIN association a ON (gp.term2_id = a.term_id)
  WHERE gp.term1_id = @t1Id
  AND a.is_not = 0
  AND gp.relationship_type_id IN (SELECT id FROM term WHERE name='part_of' OR name='is_a')
);")
dbSendQuery(con, "SET @dishin =
  ( SELECT AVG(dishin.ic)
    FROM
    (SELECT MAX(ca_ic.ic) AS ic
    FROM
    ( SELECT ca.term_id, ca.diff, -LOG(COUNT(DISTINCT a.gene_product_id)/@maxFreq) AS ic
    FROM
    ( SELECT ca.term_id,
    ABS(ca.ca_t1_number - ca.ca_t2_number) AS diff
    FROM
    (SELECT ca.ancestor AS term_id,
    COUNT(DISTINCT ca_t1_nodes.term2_id) AS ca_t1_number,
    COUNT(DISTINCT ca_t2_nodes.term2_id) AS ca_t2_number
    FROM
    ( SELECT p1.term1_id AS ancestor
    FROM graph_path p1,
    graph_path p2
    WHERE p1.term2_id = @t1Id
    AND p2.term2_id = @t2Id
    AND p1.term1_id = p2.term1_id
    AND p1.relationship_type_id IN
    (SELECT id
    FROM term
    WHERE name='part_of'
    OR name='is_a')
    AND p2.relationship_type_id IN
    (SELECT id
    FROM term
    WHERE name='part_of'
    OR name='is_a')) AS ca
    INNER JOIN graph_path ca_t1_nodes ON (ca.ancestor = ca_t1_nodes.term1_id)
    INNER JOIN graph_path ca_t2_nodes ON (ca.ancestor = ca_t2_nodes.term1_id)
    WHERE ca_t1_nodes.term2_id IN
    ( SELECT p2.term1_id AS ancestor
    FROM graph_path p2
    WHERE p2.term2_id = @t1Id)
    AND ca_t2_nodes.term2_id IN
    ( SELECT p2.term1_id AS ancestor
    FROM graph_path p2
    WHERE p2.term2_id = @t2Id)
    AND ca_t1_nodes.relationship_type_id IN
    (SELECT id
    FROM term
    WHERE name='part_of'
    OR name='is_a')
    AND ca_t2_nodes.relationship_type_id IN
    (SELECT id
    FROM term
    WHERE name='part_of'
    OR name='is_a')
    GROUP BY ca.ancestor ) AS ca ) AS ca
    INNER JOIN graph_path gp ON (ca.term_id = gp.term1_id)
    INNER JOIN association a ON (gp.term2_id = a.term_id)
    WHERE a.is_not = 0
    AND gp.relationship_type_id IN
    (SELECT id
    FROM term
    WHERE name='part_of'
    OR name='is_a')
    GROUP BY ca.term_id,
    ca.diff ) AS ca_ic
    GROUP BY ca_ic.diff) AS dishin );")
dbSendQuery(con, "SET @maxIC = ( SELECT -LOG(1/@maxFreq) );")
dbSendQuery(con, "SET @t1IC_norm = ( SELECT @t1IC/@maxIC );")
dbSendQuery(con, "SET @t2IC_norm = ( SELECT @t2IC/@maxIC );")
dbSendQuery(con, "SET @dishin_norm = ( SELECT @dishin/@maxIC );")

simRes <- dbGetQuery(con, "SELECT @dishin_norm as Sim_resnik;")
simJC  <- dbGetQuery(con, "SELECT @t1IC_norm + @t2IC_norm - 2*@dishin_norm as Dist_jc;")
simLin <- dbGetQuery(con, "SELECT (2*@dishin_norm) / (@t1IC_norm + @t2IC_norm)  as Sim_lin;")
dbDisconnect(con)

