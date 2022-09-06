var library_data = LIBRARY_JSON_DATA_PLACEHOLDER;
var table_plot = makeTable()
  .datum(library_data)
  .sortBy('Cosine', true)
  .filterCols(['Dataset', 'Status']);

d3.select('#library_table').call(table_plot);

// table_plot.on('highlight', function(data, on_off){
//   if(on_off){//if the data is highlighted
//     d3.select('#highlighted').text(
//       // data.Cosine
//     );
//   }
// });
// table_plot.on('select', function(data, on_off){
//   if(on_off){//if the data is highlighted
//     d3.select('#selected').text(
//       // data.Cosine
//     );
//   }
// });