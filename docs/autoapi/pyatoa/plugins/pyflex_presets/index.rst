:py:mod:`pyatoa.plugins.pyflex_presets`
=======================================

.. py:module:: pyatoa.plugins.pyflex_presets

.. autoapi-nested-parse::

   Pyflex requires configuration objects to set the number of available variables
   that tune the Flexwin algorithm.

   This file contains some preset maps for Pyflex that are mostly taken directly
   from the Flexwin publication Maggi et al. 2009.

   For variable descriptions see:
       https://krischer.github.io/pyflex/_modules/pyflex/config.html

   + Descriptions of a few commonly used parameters that are not self explanatory

       1. Short Term Average Long Term Average water level
           :stalta_waterlevel (float): reject windows where sta/lta waveform dips
               below this threshold value. between 0 and 1
       2. Water level rejection
           :c_0: reject if window.stalta.min() < c_0 * stalta_waterlevel
           :c_1: min_acceptable_window_length = c_1 * T_min
       3. Prominence rejection
           :c_2: reject if window.stalta.min() < c_2 * window.stalta.max()
       4. Separation height in phase separation
           :c_3a: d_stlta > c_3a * d_stalta_center * f_time
               where d_stalta = current max height above min
               and   d_stalta_center = central max height above min
               and   f_time = time decay function
           :c_3b: d_time = separation between center of window and internal maxima
               if d_time > c_3b then f_time is a time decay function, else its 1
               - if c_3b goes down,
       5. Emergent start/stops and coda wave curtailing
           :c_4a: time_decay_left = T_min * c_4a / dt
           :c_4b: time_decay_right: T_min * c_4b / dt



Module Contents
---------------

.. py:data:: pyflex_presets
   

   

