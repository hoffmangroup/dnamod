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

$(document).ready(function() {
    var idx = lunr(function() {
        this.field('Common Name');
        this.field('ChEBI Id');
        this.field('IUPAC Name');
        this.field('Synonyms');
        this.field('Chemical Formula');
        this.field('Abbreviation');
        this.field('Verified');
        this.ref('href');
    });

    var data = new Array();
    var jsonLocation = "scripts/lunr.json";
    request = new XMLHttpRequest();
    request.open('GET', jsonLocation, false);
    request.onload = function() {
        if (request.status >= 200 && request.status < 400) {
            data = JSON.parse(request.responseText);
            window.alert('JSON success')
        } else {
            window.alert('JSON error')
        }
    };
    request.send();
    
    console.log("Data Dump: ");
    for (var i in data) {
        console.log(data[i].Common_Name);
    }
       
    $('.searchbox').keyup(function() {
        var query = $(this).val();
        var result = idx.search(query);

        var resultdiv = $('#searchresults');
        if (result.length === 0) {
            resultdiv.empty()
            var append = '<li> No Results Found</li>';
            resultdiv.append(append)
        } else {
            window.alert(result.length)
            resultdiv.empty();
            resultdiv.append('<li>' + 'found' + '</li>');
        }
    });
});
              