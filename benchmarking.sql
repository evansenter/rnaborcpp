select 
  fftbor.sequence_length, 
  rnabor.time / fftbor.time, 
  fftbor.time as fftbor_avg_time, 
  fftbor.stdev_time as fftbor_stdev_time, 
  fftbor.count as fftbor_runs, 
  rnabor.time as rnabor_avg_time, 
  rnabor.stdev_time as rnabor_stdev_time, 
  rnabor.count as rnabor_runs 
from (
  select 
    algorithm, 
    sequence_length, 
    avg(time) as time, 
    std(time) as stdev_time, 
    count(*) as count 
  from
    runs 
  where 
    algorithm = 'fftbor' 
  group by 
    algorithm, 
    sequence_length
  ) as fftbor 
join (
  select 
    algorithm, 
    sequence_length, 
    avg(time) as time, 
    std(time) as stdev_time, 
    count(*) as count 
  from 
    runs 
  where 
    algorithm = 'rnabor' 
  group by 
    algorithm, 
    sequence_length
  ) as rnabor 
on fftbor.sequence_length = rnabor.sequence_length;