var library_data = LIBRARY_JSON_DATA_PLACEHOLDER;
var table_plot = makeTable()
  .datum(library_data)
  .sortBy('Cosine', false)
  .filterCols(['Dataset', 'Status']);

d3.select('#library_table').call(table_plot);

function showLibraryTable() {
  var div = document.getElementById("library_table");
  if (div.style.display === "none") {
    div.style.display = "block";
  } else {
    div.style.display = "none";
  }
}
function showParameters() {
  var div = document.getElementById("paramsDiv");
  if (div.style.display === "none") {
    div.style.display = "block";
  } else {
    div.style.display = "none";
  }
}
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