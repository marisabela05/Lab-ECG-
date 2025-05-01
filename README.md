# HRV con transformada wavelet
## Descripción 
En este laboratorio trabajamos en el análisis de la variabilidad de la frecuencia cardíaca (HRV) utilizando la transformada wavelet. El objetivo principal fue ver cómo iba cambiando la frecuencia cardíaca de la persona a lo largo del tiempo, especialmente cuando se inducían cambios que afectan la actividad simpática y parasimpática del sistema nervioso.

Para ello, adquirimos una señal ECG durante 300 segundos usando un DAQ y programamos todo en Python. Primero filtramos la señal con un filtro IIR, luego identificamos los picos R y calculamos los intervalos R-R. A partir de ahí analizamos la HRV tanto en el dominio del tiempo como en el dominio tiempo-frecuencia. En este último usamos la transformada wavelet para ver cómo se comporta la potencia espectral en diferentes bandas de frecuencia a lo largo del tiempo.

## Tener en cuenta:
1. Se usa el módulo AD8232 para la toma del ECG.
2. Frecuencia de muestreo de 2000 Hz para asegurar.
3. La duración de la adquisición es de 5 minutos.
4. Se deben tomar los tiempos para los cuales se induce el cambio en la frecuencia cardiaca.
5. Se debe instalar las librerias:
   + Nidaqmx
   + Csv.
   + Numpy.
   + Pandas.
   + Matplotlib.
   + Butter.
   + Filtfilt.
   + Pywt

6. Los electrodos deben estar bien colocados como lo indica una guía médica.
7. Evitar movimientos bruscos o cambios de posición durante la grabación de la señal.
8. Seleccionar una función wavelet adecuada para señales biológicas.

## Fundamento teórico
### Actividades del sistema nervioso

