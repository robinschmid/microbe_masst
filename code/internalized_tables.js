var library_data = LIBRARY_JSON_DATA_PLACEHOLDER;
var table_plot = makeTable("libTable")
  .datum(library_data)
  .sortBy('Cosine', false)
  .filterCols(['Dataset', 'Status']);

try {
    // table may be empty
    d3.select('#library_table').call(table_plot);
} catch(ex) {}

// add match table
var match_data = getFullMatches();
var match_table_plot = makeTable("dataTable")
  .datum(match_data)
  .sortBy('Matches', false)

d3.select('#match_table').call(match_table_plot);

function getFullMatches() {
    var data = [];
    visitAll(root, node => {
        if (hasMatches(node)) {
            var entry = {};
            entry["Name"] = node.name ?? "";
            entry["NCBI"] = node.NCBI ?? "";
            entry["Rank"] = node.Rank ?? "";
            entry["Interventions"] = node.Interventions ?? "";
            entry["Community_composition"] = node.Community_composition ?? "";
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

// add match table
var matched_datasets = getMatchedDatasets();
var dataset_table_plot = makeTable("datasetTable")
    .datum(matched_datasets)
.sortBy('Cosine', false)
// .filterCols([]);

d3.select('#dataset_table').call(dataset_table_plot);

function getMatchedDatasets() {
    var data = [];
    visitAll(root, node => {
        if (node.matches) {
            var count = node.matches.length;
            for (var i = 0; i < count; i++) {
                match = node.matches[i]
                entry = {};
                entry["Name"] = node.name ?? "";
                entry["NCBI"] = node.NCBI ?? "";
                entry["Rank"] = node.Rank ?? "";
                entry["Interventions"] = node.Interventions ?? "";
                entry["Community_composition"] = node.Community_composition ?? "";
                entry["Cosine"] = match["Cosine"] ?? "";
                entry["Matching signals"] = match["Matching Peaks"] ?? "";
                entry["Delta Mass"] = match["Delta Mass"] ?? "";
                entry["USI"] = match["USI"] ?? "";
                entry["MassIVE"] = match["USI"].split(":")[1] ?? "";
                entry["File"] = match["USI"].split(":scan:")[0] ?? "";
                data.push(entry);
            }
        }
    });
    return data;
}

function showLibraryTable() {
    toggleShowDiv("library_table");
}

function showMatchTable() {
    toggleShowDiv("match_table");
}

function showDatasetTable() {
    toggleShowDiv("dataset_table");
}

function showParameters() {
    toggleShowDiv("paramsDiv");
}

function toggleShowDiv(name) {
    var div = document.getElementById(name);
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