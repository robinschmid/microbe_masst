var library_data = LIBRARY_JSON_DATA_PLACEHOLDER;
var table_plot = makeTable("libTable")
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

// add match table
var match_data = getFullMatches();
var match_table_plot = makeTable("dataTable")
  .datum(match_data)
  .sortBy('Cosine', false)
  .filterCols([]);

d3.select('#match_table').call(match_table_plot);

function showMatchTable() {
  var div = document.getElementById("match_table");
  if (div.style.display === "none") {
    div.style.display = "block";
  } else {
    div.style.display = "none";
  }
}


function getFullMatches() {
    var data = [];
    visitAll(root, node => {
        if (hasMatches(node)) {
            var entry = {};
            entry["Name"] = node.name ?? "";
            entry["NCBI"] = node.NCBI ?? "";
            entry["Rank"] = node.Rank ?? "";
            entry["Matches"] = node.matched_size;
            entry["Samples"] = node.group_size;
            try {
                entry["Fraction"] = Number.parseFloat(node.occurrence_fraction ?? 0).toFixed(4);
            } catch(Exception) {}
            data.push(entry);
        }
    });
    return data;
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