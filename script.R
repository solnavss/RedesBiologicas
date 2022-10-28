#----------------------------------------------------------------------
# Marisol Navarro Miranda
# Cinvestav Irapuato
# Creacion: oct 2022
# Ultima modificacion: oct 2022
#----------------------------------------------------------------------

setwd("/Users/solouli/Desktop/")

library(igraph)         # Biblioteca para análisis de redes.
library(tidyverse)      # Colección de bibliotecas para manipulación de datos.
library(viridis)        # Biblioteca con paletas de colores.

# Generamos una red aleatoria de 15 nodos no dirigida.
# Para poder reproducirla, dado que tiene componentes aleatorios,
# Fijamos la semilla del generador de números pseudoaleatorios.

set.seed(25032022)
g <- sample_gnp(15, 0.2, directed = FALSE, loops = FALSE)
plot(g)

# Hacemos una nueva visualización, personalizando algunos elementos.
# https://github.com/RLadiesCuerna/meetup_2020_sep/blob/master/rcolor.pdf

set.seed(25032022)
plot(g,
     vertex.color="darkorchid",  # color de los nodos
     vertex.size=20,             # tamaño de los nodos
     edge.color="black")         # color de las aristas

# Usamos funciones de igraph para calcular:

# 1) Nodos | Vertix
V(g)

# 2) Aristas | Edges
E(g)

# 3) El grado de los nodos | Degree
degree(g)

# 4) Matriz de distancias, caminos | Distance matrix, paths
distances(g)

# 5) Componentes conexos | Connected components
is_connected(g)
count_components(g)
components(g)

# 6) Camino mas corto | Shortest path
all_shortest_paths(g, 1, to = 5)
average.path.length(g, directed=FALSE, unconnected=TRUE)

# 7) Diametro | Diameter
diameter(g)

# 8) Densidad | Density
edge_density(g)

# 9) Distribución de grado | Degree distribution
degree_distribution(g)

# 10) Centralidad,coeficiente de intermediación | Betweenness centrality
betweenness(g)

# 11) Coeficiente de clustering | Clustering coefficient
transitivity(g)

#----------------------------------------------------------------------
# OscarFontanelli
# RMB
# file:///Users/solouli/Downloads/Telegram%20Desktop/Redes_Biologicas_R.html
#----------------------------------------------------------------------

# Descargamos la lista de aristas para la red de interacciones entre proteinas.
id <- "1O5NgBBewLPHGpKDFZfNmNVu0BY-Y-dDh"
yeast <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id))
head(yeast)

# Convertimos esta lista de aristas en una red.
yeast %>%
  graph_from_data_frame -> g
# Primer vistazo a la red
g

# Nodos de la red.
V(g)

get.data.frame(g, what="vertices") %>% head()

# Aristas de la red.
E(g)

get.data.frame(g, what="edges") %>% head()

plot(g,
     vertex.size=3,
     vertex.color="blue",
     vertex.label=NA,     # No queremos ver las etiquetas de los nodos.
     edge.arrow.size=0)   # No queremos ver las puntas de las flechitas.

# Usamos la función simplify() de igraph para quitar conexiones múltiples
# y autoconexiones. Usamos la función as.unidirected() para quitar las
# direcciones de las aristas.
g <- igraph::simplify(g)
g <- as.undirected(g)
g

# Extraemos el núcleo (core) de la red formado por los nodos que tienen
# grado mayor a 3. 
core <- coreness(g,mode="all")        # Vemos a qué core pertenecen los nodos.
core <- core[core>3]                  # Extraemos los nodos del core 3.
g <- induced_subgraph(g, names(core)) # Subgrafo inducido por los nodos en el vector core.
g

V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g)
get.data.frame(g,what="vertices") %>% head() 

# Visualizamos de modo que el tamaño de cada nodo sea proporcional a 
# su grado.
plot(g,
     vertex.size=V(g)$degree,
     vertex.color="blue",
     vertex.label=NA)

# Ajustamos el tamaño de los nodos.
plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color="blue",
     vertex.label=NA)

# Lo mismo, pero el tamaño en función de la intermediación.
plot(g,
     vertex.size=log(V(g)$betweenness+1),
     vertex.color="blue",
     vertex.label=NA)

# Vemos cuáles son los nodos de mayor grado y mayor intermediación.
degree(g) %>% sort(decreasing=TRUE) %>% head(n=15)

get.data.frame(g,what="vertices") %>%
  as_tibble() %>%
  arrange(-degree) %>%
  head(n=15) %>%
  pull(name)

get.data.frame(g,what="vertices") %>%
  as_tibble() %>%
  arrange(-betweenness) %>%
  head(n=15) %>%
  pull(name)

# Búsqueda de comunidades mediante el algoritmo de Louvain.
coms <- cluster_louvain(g)
coms

coms$membership

# Agregamos a los nodos un atributo que nos indique a qué comunidad
# pertenecen.
V(g)$community <- coms$membership

# Coloreamos los nodos de acuerdo a su comunidad y probamos con diferentes
# layouts.
plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$community,
     vertex.label=NA,
     layout=layout_with_kk)

plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$community,
     vertex.label=NA,
     layout=layout_with_fr)

# Pequeño script para ver la misma red con nueve diferentes layouts.
# El último de ellos (layout_with_gem) puede tardar un poco en correr.
layouts <- c("layout_nicely","layout_with_kk","layout_in_circle",
             "layout_on_grid","layout_on_sphere","layout_randomly",
             "layout_with_fr","layout_as_star","layout_with_gem")
par(mfrow=c(3,3),mar=c(1,1,1,1))
for(i in 1:length(layouts)){
  LO <- layouts[i]
  l <- do.call(LO,list(g))
  plot(g,
       vertex.size=2*sqrt(V(g)$degree),
       vertex.color=V(g)$community,
       vertex.label=NA,
       layout=l,
       main=LO)
}

# Seleccionamos una paleta de color de la colección viridis.
# Primero visualizamos la paleta.
pal <- turbo(max(V(g)$community))
x <- seq(1,length(pal))
y <- rep(1,length(pal))
plot(x,y,cex=6,pch=19, col=pal)

# Asignamos estos colores a los nodos de la red y visualizamos.
V(g)$color <- pal[V(g)$community]
plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$color,
     vertex.label=NA,
     layout=layout_nicely)

# Primero creamos una tabla con los nodos y su comunidad.
membership.frame <- tibble(name=V(g)$name,community=V(g)$community)
# Peso de las aristas entre nodos de la misma comunidad.
edge.weight <- 50
# Lista de aristas.
edge.list <- get.edgelist(g)
# Un vector vacío para ir guardando los pesos.
weights <- c()
# Si los dos nodos de una aristas pertenecen a la comunidad,
# asignar el peso edge.weight. De lo contrario, asignar un peso de 1.
# Iterar sobre todas las aristas.
for(i in 1:nrow(edge.list)){
  com.1 <- membership.frame[membership.frame$name == edge.list[i,][1],]$community
  com.2 <- membership.frame[membership.frame$name == edge.list[i,][2],]$community
  if(com.1 == com.2){
    weight <- edge.weight
  }else{
    weight <- 1
  }
  weights <- c(weights,weight)
}
# Asignar estos pesos a las aristas.
E(g)$weight <- weights
get.data.frame(g, what="edges") %>% head(n=30)

# Visualizamos de nuevo.
set.seed(25032022)
plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$color,
     vertex.label=NA,
     layout=layout_nicely,
     edge.color="gray80")            # color de las aristas

# Creamos una tabla con las aristas. Para cada aristas nos fijamos
# en sus dos nodos. Si pertenecen a la misma comunidad, le ponemos
# a la aristas el color de los nodos. Si no, la ponemos gris.
# Iteramos sobre todas las aristas. Visualizamos.
edge.list.frame <- as.data.frame(get.edgelist(g))
E(g)$color <- map2_chr(edge.list.frame$V1, edge.list.frame$V2, ~ {
  ifelse(
    V(g)$color[V(g)[.x]] ==
      V(g)$color[V(g)[.y]],
    V(g)$color[V(g)[.x]],
    "gray70") 
})

plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$color,
     vertex.label=NA,
     layout=layout_nicely,
     edge.color=E(g)$color)

# Guardamos la red en formato graphml para su uso futuro.
write.graph(g,"D:/Bio/yeast_protein_interaction.graphml", format="graphml")
