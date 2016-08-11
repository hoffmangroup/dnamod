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
            console.log('JSON success')
        } else {
            console.log('JSON error')
        }
    };
    request.send();
    
    console.log("Data Dump: ");
    for (var i in data) {
        var doc = {
            'Common Name': data[i].CommonName,
            'ChEBI Id': data[i].ChEBIId,
            'IUPAC Name': data[i].IUPACName,
            'Synonyms': data[i].Synonyms,
            'Chemical Formula': data[i].ChemicalFormula,
            'Abbreviation': data[i].Abbreviation,
            'Verified': data[i].Verified,
            'href':data[i].CommonName + '.html'
        };
        console.log(doc)
        idx.add(doc)
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
            resultdiv.empty();
            for (var j in result) {
                console.log(result[j])
                resultdiv.append('<li>' + result[j].ref + '</li>');
            }
        }
    });
});
              