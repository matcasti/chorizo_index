
# Cargamos librerias ------------------------------------------------------

library(data.table)

# Importamos los datos ----------------------------------------------------

## Asignamos los datos a 'dt'
dt <- data.table::fread(input = "data-raw/15_a37.csv")

## Cambiamos el nombre de la columna de buceo a 'num'
names(dt)[1] <- "num"

## Asignamos datytime a POSIXct
dt[, daytime := as.POSIXct(daytime, format="%d-%m-%Y %H:%M:%S", tz = "GMT")]

## Si la profundidad es mayor o igual a 4 decimos que está buceando
dt[, is_diving := (depth >= 4)][]

## Generamos un número tipo-recorrido para cada vez que cambia el estado de buceo
dt[, dive := rleid(is_diving)][]

## Si no está buceando, asignamos `NA` al número de buceo
dt[is_diving == FALSE, dive := NA][]

## Volvemos a asignar el número de buceo para cada vez que cambia el número de
## buceo sin considerar los valores que eliminamos debido a que la profundidad
## era mayor a 4 metros
dt[!is.na(dive), dive := rleid(dive)][]

## Calculamos la diferencia del tiempo final e inicial del buceo para cada buceo
dt[!is.na(dive), diving_time := (max(daytime) - min(daytime)), dive][]

## Inspeccionamos aquellos buceos asignados
dt[!is.na(dive)]

## Eliminamos aquellos buceos que el tiempo de buceo haya sido menor a 6 segundos
dt[diving_time < 6, `:=`(dive = NA, diving_time = NA)][]

## Obtenemos el rango de tiempo de buceo
dt[!is.na(dive), range(diving_time)]

## Volvemos a asignar el número de buceo correlativo a aquellos buceos que cumplen
## con nuestros criterios (profundidad > 4m y tiempo de buceo > 6)
dt[!is.na(dive), dive := rleid(dive)][]

## Volvemos a calcular el tiempo de buceo solo para los buceos válidos
dt[!is.na(dive), diving_time := (max(daytime) - min(daytime)), dive][]

save(dt, file = "data/dt.RData")
