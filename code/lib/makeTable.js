function makeTable(id) {
    var data, sort_by, filter_cols; // Customizable variables
    var tableID = id;
    var table; // A reference to the main DataTable object

    // This is a custom event dispatcher.
    var dispatcher = d3.dispatch('highlight', 'select');

    // Main function, where the actual plotting takes place.
    function _table(targetDiv) {
        // Create and select table skeleton
        var tableSelect = targetDiv.append("table")
            .attr("class", "display compact")
            .attr("width", "100%")
            .attr("id", tableID)
            .style("visibility", "hidden"); // Hide table until style loads;

        // Set column names
        var colnames = Object.keys(data[0]);
        if (typeof filter_cols !== 'undefined') {
            // If we have filtered cols, remove them.
            colnames = colnames.filter(function (e) {
                // An index of -1 indicate an element is not in the array.
                // If the col_name can't be found in the filter_col array, retain it.
                return filter_cols.indexOf(e) < 0;
            });
        }

        // Here I initialize the table and head only.
        // I will let DataTables handle the table body.
        var headSelect = tableSelect.append("thead");
        headSelect.append("tr")
            .selectAll('td')
            .data(colnames).enter()
            .append('td')
            .html(function (d) {
                return d;
            });

        if (typeof sort_by !== 'undefined') {
            // if we have a sort_by column, format it according to datatables.
            sort_by[0] = colnames.indexOf(sort_by[0]); //colname to col idx
            sort_by = [sort_by]; //wrap it in an array
        }


        // Apply DataTable formatting: https://www.datatables.net/
        $(document).ready(function () {
            table = $('#' + tableID).DataTable({
                // Here, I am supplying DataTable with the data to fill the table.
                // This is more efficient than supplying an already contructed table.
                // Refer to http://datatables.net/manual/data#Objects for details.
                dom: 'Blfrtip',
                data: data,
                columns: colnames.map(function (e) {
                    return {
                        data: e,
                        render: function (data, type, row, meta) {
                            if (type === 'display') {
                                sdata = String(data);
                                if (sdata.startsWith("CCMSLIB"))
                                    data = '<a href="https://gnps.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID=' + data
                                        + '" target="_blank">' + data + '</a>';
                                else if (sdata.startsWith("mzspec:")) {
                                    if (sdata.includes(":scan:") || sdata.includes("accession:CCMSLIB")) {
                                        // show mirror if input usi was not empty
                                        if (inputUsi != null && inputUsi.length>0)
                                            data = '<a href="https://metabolomics-usi.gnps2.org/dashinterface/?usi1=' + inputUsi
                                                + '&usi2=' + data + '" target="_blank">' + data + '</a>';
                                        else {
                                            data = '<a href="https://metabolomics-usi.gnps2.org/dashinterface/?usi1=' +
                                                data + '" target="_blank">' + data + '</a>';
                                        }
                                    } else {
                                        data = '<a href="https://dashboard.gnps2.org/?usi=' + data
                                            + '" target="_blank">Open file</a>';
                                    }
                                } else if (sdata.startsWith("MSV0"))
                                    data = '<a href="https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=' + data
                                        + '" target="_blank">' + data + '</a>';
                            }
                            return data;
                        }
                    };
                }),
                "bLengthChange": true, // Disable page size change
                "bDeferRender": true,
                "order": sort_by,
                "lengthMenu": [2, 4, 10, 20, 50],
                "pageLength": 4
            });

            tableSelect.style("visibility", "visible");
            $('#' + tableID + ' tbody')
                .on('mouseover', 'tr', function () {
                    highlight(this, true);
                })
                .on('mouseleave', 'tr', function () {
                    highlight(this, false);
                })
                .on('click', 'tr', function () {
                    select(this);
                });
        });
    }

    /**** Helper functions to highlight and select data **************/
    function highlight(row, on_off) {
        if (typeof on_off === 'undefined') {
            // if on_off is not provided, just toggle class.
            on_off = !d3.select(row).classed('highlight');
        }
        // Set the row's class as highlighted if on==true,
        // Otherwise remove the 'highlighted' class attribute.
        // In DataTables, this is handled automatically for us.
        d3.select(row).classed('highlight', on_off);

        // Fire a highlight event, with the data and highlight status.
        dispatcher.highlight(table.rows(row).data()[0], on_off);
    }

    function select(row, on_off) {
        // Similar to highlight function.
        if (typeof on_off === 'undefined') {
            on_off = !d3.select(row).classed('selected');
        }

        d3.select(row).classed('selected', on_off);

        // Fire a select event, with the data and selected status.
        dispatcher.select(table.rows(row).data()[0], on_off);
    }

    /**** Setter / getters functions to customize the table plot *****/
    _table.datum = function (_) {
        if (!arguments.length) {
            return data;
        }
        data = _;

        return _table;
    };
    _table.filterCols = function (_) {
        if (!arguments.length) {
            return filter_cols;
        }
        filter_cols = _;

        return _table;
    };
    _table.sortBy = function (colname, ascending) {
        if (!arguments.length) {
            return sort_by;
        }

        sort_by = [];
        sort_by[0] = colname;
        sort_by[1] = ascending ? 'asc' : 'desc';

        return _table;
    };


    // This allows other objects to 'listen' to events dispatched by the _table object.
    d3.rebind(_table, dispatcher, 'on');

    // This is the return of the main function 'makeTable'
    return _table;
}