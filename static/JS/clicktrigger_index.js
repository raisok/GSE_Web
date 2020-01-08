$('#a_expdiy').click(function(){
  var gene=$( "#search-text" ).val();
  location.href= "gene=" + gene + "&clicktag=expdiy";
});
$('#a_boxplot').click(function(){
  var gene=$( "#search-text" ).val();
  location.href= "gene=" + gene + "&clicktag=boxplot";
});
$('#a_stageplot').click(function(){
  var gene=$( "#search-text" ).val();
  location.href= "gene=" + gene + "&clicktag=stageplot";
});
$('#a_survival').click(function(){
  var gene=$( "#search-text" ).val();
  location.href= "gene=" + gene + "&clicktag=survival";
});
$('#a_similar').click(function(){
  var gene=$( "#search-text" ).val();
  location.href= "gene=" + gene + "&clicktag=similar";
});
