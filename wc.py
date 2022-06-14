word_cloud_name = ""
word_cloud_des = ""
word_cloud_pap = ""
bad_words = ["infection", "Infection", "human", "Human", "pathogen", "Pathogen",
             "Bacterium", "bacterium", "Disease", "disease", "Clinical", "clinical"]
for bact_tmp in bacts:
    name = bact_tmp.name
    des = bact_tmp.description
    pap = ""
    for paper in bact_tmp.papers:
        pap += paper.title
    for bw in bad_words:
        name = name.replace(bw, "")
        des = des.replace(bw, "")
        pap = pap.replace(bw, "")
    word_cloud_name += (name + " ")
    word_cloud_des += (des + " ")
    word_cloud_pap += (pap + " ")

wordcloud = WordCloud(max_font_size=40, background_color="white", contour_color='#5e81ac',
                      contour_width=0.1, scale=2).generate(remove_stopwords(word_cloud_name))
plt.figure()
plt.imshow(wordcloud, interpolation="bilinear")
plt.axis("off")
plt.savefig("plot_print/name_wc.png", dpi=500)
plt.close()

wordcloud = WordCloud(max_font_size=40, background_color="white", contour_color='#5e81ac',
                      contour_width=0.1, scale=2).generate(remove_stopwords(word_cloud_des))
plt.figure()
plt.imshow(wordcloud, interpolation="bilinear")
plt.axis("off")
plt.savefig("plot_print/des_wc.png", dpi=500)

wordcloud = WordCloud(max_font_size=40, background_color="white", contour_color='#5e81ac',
                      contour_width=0.1, scale=2).generate(remove_stopwords(word_cloud_pap))
plt.figure()
plt.imshow(wordcloud, interpolation="bilinear")
plt.axis("off")
plt.savefig("plot_print/pap_wc.png", dpi=500)
