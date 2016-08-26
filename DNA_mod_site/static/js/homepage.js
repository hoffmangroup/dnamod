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
    $("#tabs").tabs({
        active:1
    });
    $('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
        $('.searchbox').val('');
        $('.searchbox').focus().blur();
        $('#searchresults').empty();
    });
});

$(document).ready(function() {
    var idx = elasticlunr(function() {
        this.addField('Common Name');
        this.addField('ChEBI Id');
        this.addField('IUPAC Name');
        this.addField('Synonyms');
        this.addField('Chemical Formula');
        this.addField('Abbreviation');
        this.addField('Verified');
        this.addField('Symbol');
        this.setRef('href');
        this.saveDocument(false);
    });

    var data = new Array();
    var jsonLocation = "js/lunr.json";
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
    
    var store = new Array();
    for (var i in data) {
        var doc = {
            'Common Name': data[i].CommonName,
            'ChEBI Id': data[i].ChEBIId,
            'IUPAC Name': data[i].IUPACName,
            'Synonyms': data[i].Synonyms,
            'Chemical Formula': data[i].ChemicalFormula,
            'Abbreviation': data[i].Abbreviation,
            'Verified': data[i].Verified,
            'Symbol': data[i].Symbol,
            'href':data[i].Refname + '.html'
        };

        store[doc.href] = {
            'Common Name': data[i].CommonName,
            'ChEBI Id': data[i].ChEBIId,
            'IUPAC Name': data[i].IUPACName,
            'Synonyms': data[i].Synonyms,
            'Chemical Formula': data[i].ChemicalFormula,
            'Abbreviation': data[i].Abbreviation,
            'Verified': data[i].Verified,
            'Symbol': data[i].Symbol
        };
        
        idx.addDoc(doc)
    }
       
    $('.searchbox').keyup(function() {
        var query = $(this).val();
        
        if (query.length <= 2) {
            var searchfields = ['Abbreviation', 'Symbol'];
        } else {
            var searchfields = ['Common Name', 'Abbreviation', 'Synonyms', 'Symbol', 'ChEBI Id', 'IUPAC Name', 'Chemical Formula'];
        }
        
        var resultdiv = $('#searchresults');
        var displayedResults = new Array;
        
        resultdiv.empty();
        for (var field in searchfields) {
            var category = searchfields[field];
            var numboost = 1;
            
            if (query.length <= 2) {
                numboost = 10;
            }
            
            var configuration = '{"fields": {"'+category+'": {"boost": "'+numboost+'"}}, "expand": true}';
            configuration = JSON.parse(configuration);
            
            var result = idx.search(query, configuration);
            
            var uniqueCount = 0;
            for (i in result) {
                inArray = $.inArray(store[result[i].ref]['Common Name'], displayedResults);
                if (inArray === -1) {
                    uniqueCount = uniqueCount + 1;
                }
            }
            
            if (uniqueCount > 0) {
                resultdiv.append('<h3> Query matches: ' + category + ' </h3>');
            }
            
            for (var j in result) {
                inArray = $.inArray(store[result[j].ref]['Common Name'], displayedResults);
                if (inArray === -1) {
                    displayedResults.push(store[result[j].ref]['Common Name']);
                    if (store[result[j].ref]['Verified'] === 1) {
                        resultdiv.append('<li class="list-group-item"> <a href= "' + result[j].ref + '" style="color:#4DAC26">' + store[result[j].ref]['Common Name'] + '</a> </li>');
                    } else {
                        resultdiv.append('<li class="list-group-item"> <a href= "' + result[j].ref + '" style="color:#D01C8B">' + store[result[j].ref]['Common Name'] + '</a> </li>');
                    }
                }
            }
        }
    });
    
});