Taller: Análisis de redes biológicas
================
Marisol Navarro-Miranda
28/10/22

**Análisis de redes biológicas**

# Temas

1 Descripción de la red

* V (nodos)
* E (aristas)
* Grado
* Caminos
* Componentes conexas
* Camino más corto
* Diámetro
* Densidad
* Distribución de grado
* Centrality (Degree y Betweenness)
* Coeficiente de clustering

2 Estudio de caso

# 1. Descripción de la red

**R** 

``` r
#----------------------------------------------------------------------
# Marisol Navarro Miranda
# Cinvestav Irapuato
# Creacion: oct 2022
# Ultima modificacion: oct 2022
#----------------------------------------------------------------------

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
```
# 2. Estudio de caso

**R** 

``` r
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

# Visualizamos de nuevo. Pero esta vez salvamos.

png("network.png", width = 300*10, height = 300*8,
    res = 300, units = "px")

set.seed(25032022)
plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$color,
     vertex.label=NA,
     layout=layout_nicely,
     edge.color="gray80")            # color de las aristas

dev.off()

# Guardamos la red en formato graphml para su uso futuro.

write.graph(g,"/Users/solouli/Desktop/yeast_protein_interaction.graphml", format="graphml")
```

# Referencias

* Recursos

https://igraph.org/r/html/latest/
http://networksciencebook.com/chapter/1

https://kateto.net/2016/05/network-datasets/
https://kateto.net/wp-content/uploads/2016/01/NetSciX_2016_Workshop.pdf
https://github.com/elaragon/R-igraph-Network-Workshop/blob/master/NetSciX%202016%20Workshop.R

* Bases de Datos

https://snap.stanford.edu/data/
https://networkrepository.com/
https://icon.colorado.edu/#!/networks

* Crea un dataset

https://www.genecards.org/
https://string-db.org/cgi/input?sessionId=bvevnhF1MzII&input_page_show_search=on



