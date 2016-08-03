$(document).ready(function() {
    $("#accordion").accordion({
    heightStyle: "content",
    collapsible: "true",
    active: false,
    });
    $("#accordion2").accordion({
    heightStyle: "content",
    collapsible: "true",
    active: false,
    });
    $("#tabs").tabs();
});

$(document).ready(function () {
    var index, store;
    $.getJSON('/lunr.json', function(response) {
        index = lunr.Index.load(response.index);
        store = response.store;
        
        $('.searchbox').on('keyup', function() {
            var query = $(this).val();
            var result = index.search(query);
            var resultdiv = $('#searchresults');
            
            if (result.length === 0) {
                resultdiv.hide();
            } else {
                resultdiv.empty();
                for (var item in result) {
                    var ref = result[item];
                    var searchitem = 'result';
                    console.log(searchitem)
                }
                resultdiv.append(searchitem);
            }
        });
    });
});
                    