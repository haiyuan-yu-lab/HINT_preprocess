NEW_DIR='/home/yl986/data/HINT/update_2026/outputs/202605_1'
PREV_DIR='/home/yl986/data/HINT/update_2024/outputs/202410_1'
# PREV_DIR='/home/yl986/data/HINT/previous_versions/2021-08'
python -u ../5-plot_venn.py \
    --hint_new_dir $NEW_DIR \
    --hint_prev_dir $PREV_DIR \
    --timetag_new '202605_rev' \
    --timetag_old '202410_rev'


# NEW_DIR='/home/yl986/data/HINT/update_2024/outputs/202410_1'
# # PREV_DIR='/home/yl986/data/HINT/update_2024/outputs/202410'
# PREV_DIR='/home/yl986/data/HINT/previous_versions/2021-08'
# python -u ../5-plot_venn.py \
#     --hint_new_dir $NEW_DIR \
#     --hint_prev_dir $PREV_DIR \
#     --timetag_new '202410_rev' \
#     --timetag_old '202108'